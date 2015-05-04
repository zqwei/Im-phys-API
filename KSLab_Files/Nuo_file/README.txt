What is included:

- Processed data files are in the folder "datafiles".  

	Each file "data_structure_XXX.mat" contains the data from one expertimental session 
	(e.g. data_structure_ANM210861_20130701.mat contains data from animal 210861 collected on 2013/07/01)
	
	Each file "meta_data_XXX.mat" contains the meta data information for each session. 
	There is one "meta_data_XXX.mat" file for each "data_structure_XXX.mat" file.

	The description of the data structures are in .\documentations\data description_NL_revision20140905_rev2




- A description of the processed data structure is in the folder "documentations"

	A couple of simple demo scripts are included that performs some basic analyses of the data
	see "data description_NL_revision20140905_rev2" for description



- A collection of analyses scripts that reproduces the figures in "Li, Chen, Guo, Gerfen, Svoboda (2015)" is included in the folder "analyses_scripts"
	
	To run a script, (e.g. "analysis_ALM_population_selectivity"):
	1) open MATLAB
	2) change current folder to the ".\analyses_scripts\" folder
	3) type run('analysis_ALM_population_selectivity')


