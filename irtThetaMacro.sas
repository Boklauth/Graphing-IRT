/*
irtTheta creates a graph for dichotomous IRT test information 
and standard error of estimate (SEE) by person location (theta);
and a graph for log likelihood function by score pattern. 
Copyright (C) 2019  Bo Klauth

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
********************************;
/*
ds= enter your dataset name
model= enter your model the models available are onep, twop, threep, and fourp.
start_var= the name of your first item
end_var=the name of your last item
nitems= the total number of items
method=ML
ods_graphics= off or on. "on" will produce item characteristic curves;
*/


*Exmaple;
/*%Theta(ds=SampleData2, model=twop, start_var=item1, end_var=item8, nitems=8, 
	method=ML, ods_graphics=off);*/
* Note that the method is always ML.;


%macro irtTheta (ds, model, start_var, end_var, nitems, method, ods_graphics);

* data step to get the composite score;
data temp;
set &ds;
total = sum (of &start_var--&end_var);
new_id=_n_;
run;

%let dataset=temp;

PROC SQL; 
select nobs INTO: nobs 
from sashelp.vtable
where libname="WORK" and memname="TEMP"
;
quit;

* Using proc irt to get the output;
ods graphics &ods_graphics;
proc irt data=&dataset itemstat itemfit 
		 scoremethod=&method out=theta outmodel=model_output
		 plots=(ICC(unpack) IIC(unpack) TIC);
	var &start_var--&end_var;
	model &start_var--&end_var
			/resfunc=&model;
	title1 "&model IRT";
	run;
	ods graphics off;

* Transforming data to get ceiling(d), guessing(c), intercept, and slope;
data want(keep=_SUBTYP_ _NAME_ _VAR1_ _EST_ _STDERR_);
set model_output;
run;
data ceiling (drop = _name_ rename=(_est_=d));
set want;
if _SUBTYP_="CEILING" & _VAR1_^="";
run; 

data guessing (drop = _name_ rename=(_est_=c));
set want;
if _SUBTYP_="GUESSING" & _VAR1_^="";
run; 

data intercept (drop = _name_ rename=(_est_=intercept));
set want;
if _SUBTYP_="INTERCEPT";
run; 

data slope (drop=_name_ rename=(_est_=slope));
set want;
if _SUBTYP_="SLOPE" & _VAR1_^="";
run;  

* Using the intercept and slope to find the item difficulty/item location;
data item_diff;
merge ceiling guessing intercept slope;
item_diff=-intercept/slope;
run;
data item_diff (drop=_subtyp_ _STDERR_);
set item_diff;
run;

* Printing all parameter estimates;
proc print data=item_diff;
label _var1_ ="Item" item_diff="p" slope="Slope" intercept="Intercept";
title "Guessing, Item Difficult and Item Discrimination Estimates for &model Model";
footnote '';
run;


* Transposing data;
proc transpose data=ceiling out=tceiling prefix=d;
var d;
run;

proc transpose data=guessing out=tguessing prefix=c;
var c;
run;

proc transpose data=slope out=tslope prefix=slope;
var slope;
run;

proc transpose data=item_diff out=titem_diff prefix=p;
var item_diff;
run;

data slope_p;
merge tceiling tguessing tslope titem_diff;
run;

*--------------------------------;
*This is to calculate total test information and SEE;
/* 
Item information only concerns with the probability of getting each item correct
from -4 to +4. It does not involve score pattern. 
And the thus, thus n = (range/increment of i) + 1= 801
*/
data Iteminfo_SEE (rename=(i=Theta));
set slope_p;
do i = -4 to 4 by 0.01;
	array d[&nitems] d1-d&nitems;
	array c[&nitems] c1-c&nitems;
	array alpha[&nitems]slope1-slope&nitems;
	array item_diff[&nitems] p1-p&nitems;
	array pj[&nitems] pj1-pj&nitems;
	array iteminfo[&nitems]iteminfo1-iteminfo&nitems;
	do j=1 to &nitems;
	pj[j]=c[j]+(d[j]-c[j])*1/(1+exp(-alpha[j]*(i-item_diff[j])));

	* item info and SEE for onep to twop models;
	%if &model=onep or &model=twop %then %do;
	iteminfo[j]=(alpha[j]**2)*(pj[j])*(1-pj[j]);
	%end;

	*item info and see for threep model;
	%if &model=threep %then %do;
	iteminfo[j]=(alpha[j]**2)*((pj[j]-c[j])**2)*(1-pj[j])/((1-c[j])**2)*(pj[j]);
	%end;

	*item info and see for fourp model;
	%if &model=fourp %then %do;
	iteminfo[j]=((alpha[j]**2)*((pj[j]-c[j])**2)*(d[j]-pj[j])**2)/
		(((d[j]-c[j])**2)*(pj[j])*(1-pj[j]));
	%end;

	iteminfototal=sum( of iteminfo1-iteminfo&nitems);
	SEE=sqrt(1/iteminfototal);
	drop j;
	end;
output;
end;

*--------------------------------;

*build Graph SEE and test info;
axis1 label=(a=90 h=2 'Test Information');
axis2 label=(h=2 'Theta');
Axis3 label=(a=90 h=2 'SEE');

symbol1 w = 2 v=none color=blue interpol=join;
symbol2 w = 2 v=none color=red interpol=join;
legend1 label =(height=1.5 'Test Information') value =('');
legend2 label=(height=1.5 'SEE') value=('');

* graphing item information by SEE;
proc gplot data=Iteminfo_SEE;
title "Test Information and SEE Plot for &model Model";
plot1 iteminfototal*Theta/hminor=1 vaxis=axis1 haxis=axis2 legend=legend1;
plot2 SEE*Theta /vaxis=axis3 legend=legend2;
run;
quit;


* preparing data steps for calcualting irt score;
proc means data=theta mean std noprint;
var total;
output out= z;
run;


data z(keep=_stat_ total);
set z;
if _stat_="N" or _stat_="MIN" or _stat_="MAX" then delete;
run;

proc transpose data=z out=z_tr;
id _stat_;
run;


data z_tr_obs (drop=_name_ i);
set z_tr;
do i = 1 to &nobs;
output;
end;
run;

data IRT_score;
merge theta z_tr_obs;
run;
* calculating IRT score;
data IRTscore;
set IRT_score;
IRT_score=MEAN + _Factor1 * STD;
new_id = _n_;
run;

*Obtaining slopes and item difficulties and prepare them for all the observations;

data itemdiff_slope (keep=d1-d&nitems c1-c&nitems slope1-slope&nitems p1-p&nitems); 
set iteminfo_see;
new_id=_n_;
if new_id=1;
do i=1 to &nobs;
output;
end;
run;
* Merge tables for items and slopes and item difficulties togeter;
data wantxxx; 
merge &dataset itemdiff_slope;
run;
* CALCULATING PERSON PROBABILITY OF EACH ITEM;
/* 
Calculating probability of each item for every observation from person locations
ranged from -i to +i;
* So each probability (pj1 to pjn) is the probability of getting an item correct;
* This is to prepare for calculating the pattern probability 
or the log-likelihood function.
*/
data wantxx (rename=(i=Theta));
set wantxxx;
do i = -4 to 4 by 0.01;
array dj[&nitems] d1-d&nitems;
array cj[&nitems] c1-c&nitems;
array slope[&nitems]slope1-slope&nitems;/*alpha*/
array p[&nitems]p1-p&nitems;
array pj[&nitems]pj1-pj&nitems;
do j = 1 to &nitems;
/*
adding dj cj in the pj calculation to make a generic calculation
for all four models
*/
pj[j]=cj[j]+(dj[j]-cj[j])*1/(1+exp(-slope[j]*(i-p[j]))); 
drop j;
end;
output;
end;
run;




* CALCULATING PROBABILITY OF A SCORE PATTERN OR LLF OF A SCORE PATTERN;
data wantx;
set wantxx;
/*array for test items*/
array items[&nitems]&start_var--&end_var;

/* array for probability of getting an item correct*/
array pj[&nitems]pj1-pj&nitems; 

/* array for probability of item according the actual response*/
array phi[&nitems] phi1-phi&nitems;

/*arrya for log likelihood of each item based on the actual response*/
array Ln[&nitems] Ln1-Ln&nitems; 
do K= 1 to &nitems;

if items[K]=1 then phi[k]=pj[K];
if items[K]=0 then phi[k]=1-pj[K];

/* Total probabilty for an examinee = multiplication of probability of each item, 
which is equivalent to using geomean*/
total_phi=geomean( of phi1-phi&nitems)**n(of phi1-phi&nitems); 
Ln[K]=log(phi[K]);
Total_Ln=sum (of Ln1-Ln&nitems);
/*checking phi (probability) and log calculation*/
Exp=exp(total_Ln);
drop k;
end;
run;

*Creating score pattern;
data wantx1;
set wantx;
length score_pattern $&nitems;
array score &start_var--&end_var;
do over score;
if not missing(score) then score_pattern=cats(score_pattern, score);
else score_pattern=catt(score_pattern, score);
pattern_length=length(score_pattern); /*checking the length of the pattern*/
end;
run;

data AcrossContinuum_Parameters;
set wantx1;
run;

* Creating frequency for score pattern;
*-----------;

*checking frequency of total score and response pattern;
* Converting simulatign frequency to dataset frequency;
proc freq data=Acrosscontinuum_parameters noprint;
Tables total /out =simTotalScore;
run;

proc freq data=Acrosscontinuum_parameters noprint;
Tables score_pattern /out =simScorePattern;
run;

data OriTotalScore;
set simTotalScore;
Frequency=COUNT/((8/0.01)+1);
label total="Total Score";
run;

data OriScorePattern;
set simScorePattern;
Frequency=COUNT/((8/0.01)+1);
label score_pattern="Score Pattern";
run;

proc print data=oriTotalScore;
var total Frequency percent;
label total="Total Score";
title "Frequency and Percentage of Total Composite Score";
run;

proc print data=oriScorePattern;
var score_pattern Frequency percent;
title "Frequency and Percentage of Score Pattern";
run;


*creating data for log likelihood function graphs;
proc sql; 
create table graph_loginfo as 
select distinct(score_pattern),total, theta, total_ln
from wantx1;
quit;

proc sort data=graph_loginfo;
by total;
run;

*build log likelihood function (LLF) and score pattern;

axis1 label=(a=90 h=2 'Log Likelihood');
axis2 label=(h=2 'Theta');

symbol1 color=black interpol=join r=0;
symbol2 color=blue interpol=join r=0;
symbol3 color=green interpol=join r=0;
symbol4 color=brown interpol=join r=0;
symbol5 color=gray interpol=join r=0;
symbol6 color=green interpol=join r=0;
symbol7 color=red interpol=join r=0;
symbol8 color=blue interpol=join r=0;
symbol9 color=gray interpol=join r=0;
symbol10 color=black interpol=join r=0;
symbol11 color=red interpol=join r=0;
symbol12 color=red interpol=join r=5;
symbol13 color=orange interpol=join r=3;
symbol14 color=pink interpol=join r=2;
symbol15 color=red interpol=join r=3;
symbol16 color=green interpol=join r=13;
symbol17 color=blue interpol=join r=13;
symbol18 color=blue interpol=join r=5;
symbol19 color=red interpol=join r=13;
symbol20 color=black interpol=join r=5;
symbol21 color=blue interpol=join r=10;
symbol22 color=black interpol=join r=1000;



proc gplot data=graph_loginfo;
plot Total_Ln*Theta=score_pattern/ vaxis=axis1 haxis= -4 to 4;
by total;
/*where total^=0;*/
title "Log Likelihood Function for Total Score Using &model Model";
run;
quit;

*-----------;

* Sort the data set by id;
proc sort data=wantx1;
by new_id;
run;


* Finding max of log likelihood;
proc sort data=wantx1;
by new_id Total_Ln;
run;

* Select only records with maximum likelihood;
proc sql;
create table max_ln as
select *,max(total_ln) as max_tln
from wantx1
group by new_id
having max(total_ln)=total_ln;
quit;

* working on here;
proc sort data=irt_score;
by new_id;
run;

data _Factor1(keep=_Factor1);
set irt_score;
run;

Data Dataset_Parameters;
merge max_ln _Factor1;
run;

Data Comparing_theta(keep=new_id total theta _Factor1 total_ln total_phi exp);
set Dataset_Parameters;
run;

* Printing 45 observations to compare the estimates of the macro and SAS;
proc print data=Comparing_theta (obs=45);
title "Comparing Theta for &model  (printed 45 observations)";
var new_id total theta _factor1;
run;

* deleting unused datasets;
proc datasets lib=work nolist;  
		delete  want wantx wantx1;
/*		delete graph_loginfo;*/
		delete wantxx;
		delete wantxxx;
		delete z z_tr z_tr_obs;
		delete guessing intercept;
		delete item_diff;
		delete slope;
		delete slope_p;
		delete itemdiff_slope;
		delete tguessing titem_diff;
		delete tslope ceiling tceiling;
		delete max_ln Temp _Factor1 Theta irt_score;
quit; 
run;

%Mend IRTTheta;

