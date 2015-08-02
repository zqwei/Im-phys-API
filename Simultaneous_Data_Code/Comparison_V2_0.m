%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  generate data from shuffled dataset code (based on
%     Comparison_V1_0 && Comparison_V1_0_1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V2_0

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
DataListShuffle        = DataSetList;

% Rename the DataSetList Names as the simultaneously recorded datasets
DataSetList(1).name    = 'Simultaneous_Spikes';
DataSetList(2).name    = 'Simultaneous_Ca_Fast_Short_Delay';
DataSetList(3).name    = 'Simultaneous_Ca_Slow_Short_Delay';
DataSetList(4).name    = 'Simultaneous_Ca_Fast_Long_Delay';
DataSetList(5).name    = 'Simultaneous_Ca_Slow_Long_Delay';
DataSetList(6).name    = 'Simultaneous_Ca_Slow_Short_Delay_Virus';
% this should be identical to Spikes analysis results after normalization
DataSetList(7).name    = 'Simultaneous_DFFSpikes'; 

minUnitsSession        = 3;
minRate                = 5;
perMinRate             = 0.4; % percentage of min Rate in an entire trial
% In the following analysis of simultaneous recorded datasets, we only use
% a subset of strong selective neurons (a trick to decrease the number of
% neurons in each session)
ROCThres               = 0.50; 

%%% 
% Case 1:
% Taking the fact that there could be different numbers of trials in for 
% units in the same ephys session

%%%
% Plot trials

% nData                  = 1;
% load([TempDatDir DataListShuffle(nData).name '.mat']);
% plotSimultaneousSpikeData(nDataSet(DataSetList(nData).ActiveNeuronIndex),...
%     DataSetList(nData).params, minRate, perMinRate, ROCThres, minUnitsSession)

for nData                  = 1;
    load([TempDatDir DataListShuffle(nData).name '.mat']);
    [nDataSet, DataSetList(nData).kickOutIndex] = ...
                            getSimultaneousSpikeData(nDataSet, ...
                            DataSetList(nData),...
                            minRate, perMinRate, ROCThres, ...
                            minUnitsSession); 
    save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');
end

%%% 
% Case 2:
% Taking the fact that there is the same number of trials in for units in
% the indentical Ca++ imaging session

for nData                  = 2:length(DataSetList)-1
    load([TempDatDir DataListShuffle(nData).name '.mat']);
    [nDataSet, DataSetList(nData).kickOutIndex] = ...
                            getSimultaneousCaimagingData(nDataSet, ...
                            DataSetList(nData), ...
                            ROCThres, minUnitsSession);  
    save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');    
end

save([TempDatDir 'DataListSimultaneous.mat'], 'DataSetList');



