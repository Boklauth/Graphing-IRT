*DEMONSTRATING ON HOW TO USE THE IRTTheata Macro;
* -----------------------------;


 
* HOW TO USE THE MACRO;
*Include macro by specifying its location on your computer;

%include "C:\Users\Dell\Google Drive\PROGRAMMING\SAS Macro\Find probability of item patterns\IRTTheta_Macro\IRTMacroPackage\irtThetaMacro.sas";

* Enter dataset and its information for the macro;

/*Simply just fill out the variables below.*/
	/*%IRTTheta(ds=, model=, start_var=, end_var=, nitems=, method=, ods_graphics=off);*/
		/**/
			/*ds= enter your dataset name*/
			/*model= enter your model; the models available are onep, twop, threep, and fourp.*/
			/*start_var= the name of your first test item in your dataset.*/
			/*end_var=the name of your last test item in your dataset.*/
			/*nitems= the total number of test items*/
			/*method=ML (Marginal Maximum Likelihood)*/
			/*ods_graphics= off or on; "on" will produce item characterisic curves using proc irt;*/

* Be sure that the dataset is accessible by the macro;
* EXAMPLE WITH THE PROVIDED DATASET;

%irtTheta(ds=sampleData2, model=threep, start_var=item1, end_var=item8, nitems=8, method=ML, ods_graphics=off);



