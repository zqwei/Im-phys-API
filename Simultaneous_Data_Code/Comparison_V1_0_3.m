%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  generate data from shuffled dataset code (based on
%     Comparison_V1_0 && Comparison_V1_0_1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_0_3

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat'], 'DataSetList', 'ephysCellIndex', 'fileToAnalysis');
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



