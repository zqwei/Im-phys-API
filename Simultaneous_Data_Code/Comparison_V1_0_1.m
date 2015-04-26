% 
% Comparison based on single unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0
%
% Comparison list
%
% 1.  Raw activity of all neurons sorted in different ways
% 2.  Shape of activity distribution across trial periods (pre-sample,
%     sample, delay, response)
% 3.  Selectivity over time
% 4.  Single neuron choice probability index distribution (ROC), across
%     trial periods
% 5.  Per session decodability over time
% 6.  Collected population decision decodability over time
% 7.  Saturation of decodeability with number of neurons
% 8.  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched
% 9.  Measures of trial-to-trial variability, by trial period
% 10. Fraction of variance captured by PCA, per session, by trial period.
%
% -------------------------------------------------------------------------
% version 1.1
% 
% + Save the dataset while before rerun all the code. 
% Simultaneously recording data
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Shape of activity distribution across trial periods (pre-sample,
%     sample, delay, response)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Selectivity over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Single neuron choice probability index distribution (ROC), across
%     trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.  Per session decodability over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.  Collected population decision decodability over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_6
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.  Saturation of decodeability with number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_7
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_8
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9.  Measures of trial-to-trial variability, by trial period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_9
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Fraction of variance captured by PCA, per session, by trial period.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_10
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat'], 'DataSetList', 'ephysCellIndex', 'fileToAnalysis', '');
DataListShuffle        = DataSetList;
fileList               = {SpikeFileList; CaImagingShortDelayFastFileList; CaImagingShortDelaySlowFileList;...
                            CaImagingLongDelayFastFileList; CaImagingLongDelaySlowFileList; CaImagingShortDelaySlowVirusFileList};
clear DataSetList;
DataSetList(1).name    = 'Spikes';
DataSetList(2).name    = 'Ca_Fast_Short_Delay';
DataSetList(3).name    = 'Ca_Slow_Short_Delay';
DataSetList(4).name    = 'Ca_Fast_Long_Delay';
DataSetList(5).name    = 'Ca_Slow_Long_Delay';
DataSetList(6).name    = 'Ca_Slow_Short_Delay_Virus';

minRate                = 5;
perMinRate             = 0.4;
% ROCThres               = 0.60;
minUnitsSession        = 3;


nData                  = 1;
load([TempDatDir DataListShuffle(nData).name '.mat']);
[nDataSet3D, nDataSet] = getSimultaneousSpikeData(nDataSet(ephysCellIndex{nData}), DataListShuffle(nData).params, minRate, perMinRate, ROCThres, minUnitsSession);  
DataSetList(nData).params  = DataListShuffle(nData).params; 
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet', 'nDataSet3D');

for nData              = 2:length(DataSetList)
    load([TempDatDir DataListShuffle(nData).name '.mat']);
    [nDataSet3D, nDataSet] = getSimultaneousCaimagingData(nDataSet(ephysCellIndex{nData}), DataListShuffle(nData).params, ROCThres, minUnitsSession);  
    DataSetList(nData).params  = DataListShuffle(nData).params; 
    save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet', 'nDataSet3D');    
end

save([TempDatDir 'DataList.mat'], 'DataSetList', 'fileToAnalysis');

