% Data Analysis All

disp('Slicing the analysis cubes in each dataset')

load ('TempDat/DataListShuffle.mat');
nData                    = 1;
load(['TempDat/' DataSetList(nData).name '.mat']);
ephysCellIndex{nData}        = [DataSetList(nData).cellinfo(:).depth]  > 120 &...
                           [DataSetList(nData).cellinfo(:).depth]  < 750 & ...
                           [DataSetList(nData).cellinfo(:).ML_axis]  > 1100 &...
                           [DataSetList(nData).cellinfo(:).ML_axis]  < 1900;                        
nData           = 2;
load(['TempDat/' DataSetList(nData).name '.mat']);
ephysCellIndex{nData}  = [DataSetList(nData).cellinfo(:).depth]  > 120 &...
                       [DataSetList(nData).cellinfo(:).depth]  < 750 & ...
                       [DataSetList(nData).cellinfo(:).AP_axis]  > 2450 &...
                       [DataSetList(nData).cellinfo(:).AP_axis]  < 2800 &...
                       [DataSetList(nData).cellinfo(:).ML_axis]  > 1100 &...
                       [DataSetList(nData).cellinfo(:).ML_axis]  < 1900; 
nData           = 3;
load(['TempDat/' DataSetList(nData).name '.mat']);
ephysCellIndex{nData}  = [DataSetList(nData).cellinfo(:).depth]  > 120 &...
                       [DataSetList(nData).cellinfo(:).depth]  < 750 & ...
                       [DataSetList(nData).cellinfo(:).AP_axis]  > 2100 &...
                       [DataSetList(nData).cellinfo(:).AP_axis]  < 2800 &...
                       [DataSetList(nData).cellinfo(:).ML_axis]  > 1100 &...
                       [DataSetList(nData).cellinfo(:).ML_axis]  < 1900; 
nData           = 4;
load(['TempDat/' DataSetList(nData).name '.mat']);
ephysCellIndex{nData}  = [DataSetList(nData).cellinfo(:).depth]  > 120 &...
                       [DataSetList(nData).cellinfo(:).depth]  < 750 & ...
                       [DataSetList(nData).cellinfo(:).AP_axis]  > 2100 &...
                       [DataSetList(nData).cellinfo(:).AP_axis]  < 2800 &...
                       [DataSetList(nData).cellinfo(:).ML_axis]  > 1100 &...
                       [DataSetList(nData).cellinfo(:).ML_axis]  < 1900;                        

nData           = 5;
load(['TempDat/' DataSetList(nData).name '.mat']);
ephysCellIndex{nData}  = [DataSetList(nData).cellinfo(:).depth]  > 120 &...
                       [DataSetList(nData).cellinfo(:).depth]  < 750 & ...
                       [DataSetList(nData).cellinfo(:).AP_axis]  > 2700 &...
                       [DataSetList(nData).cellinfo(:).AP_axis]  < 3000 &...
                       [DataSetList(nData).cellinfo(:).ML_axis]  > 1100 &...
                       [DataSetList(nData).cellinfo(:).ML_axis]  < 1900; 
nData           = 6;
load(['TempDat/' DataSetList(nData).name '.mat']);
ephysCellIndex{nData}  = [DataSetList(nData).cellinfo(:).depth]  > 120 &...
                       [DataSetList(nData).cellinfo(:).depth]  < 750 & ...
                       [DataSetList(nData).cellinfo(:).AP_axis]  > 2100 &...
                       [DataSetList(nData).cellinfo(:).AP_axis]  < 2800 &...
                       [DataSetList(nData).cellinfo(:).ML_axis]  > 1100 &...
                       [DataSetList(nData).cellinfo(:).ML_axis]  < 1900;

ROCThres                      = 0.75;
                   
fileToAnalysis                = [1 3 4 5 6]; % 1:length(DataSetList);                   
save ('TempDat/DataListShuffle.mat', 'ephysCellIndex', 'fileToAnalysis', '-append');

disp('Generating simultaneous activity dataset')
run Simultaneous_Data_Code/Comparison_V1_0_1.m

disp('Generating summary of selectivity onset time for each dataset')
run Shuffled_Data_Code/Comparison_V1_0_7.m
run Simultaneous_Data_Code/Comparison_V1_0_2.m

numRandPickUnits              = 30;

disp('Generating figures for neuronal activities')
run Shuffled_Data_Code/Comparison_V1_1.m
run Simultaneous_Data_Code/Comparison_V1_1_1.m

disp('Generating figures for neuronal activities distribution')
run Shuffled_Data_Code/Comparison_V1_2.m

disp('Generating figures for ROC')
run Shuffled_Data_Code/Comparison_V1_4.m

disp('Generating figures for selectivity of trial type')
% run Simultaneous_Data_Code/Comparison_V1_5.m
% run Simultaneous_Data_Code/Comparison_V1_6.m
run Simultaneous_Data_Code/Comparison_V1_6_1.m

% disp('Generating figures for selectivity of trial type after KO')
% run Simultaneous_Data_Code/Comparison_V1_7.m
run Simultaneous_Data_Code/Comparison_V1_7_1.m
% 
disp('Generating figures for temporal precision')
% run Simultaneous_Data_Code/Comparison_V1_8.m
run Simultaneous_Data_Code/Comparison_V1_8_1.m
% run Simultaneous_Data_Code/Comparison_V1_8_2.m


% disp('Generating figures for LDA-PCA')
run Simultaneous_Data_Code/Comparison_V1_10.m
% 
% disp('Generating figures for dPCA')
run Simultaneous_Data_Code/Comparison_V1_11.m