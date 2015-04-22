% 
% Comparison based on single unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0.4
% 
% 
% 
% Check all neurons stastical properity comparing to ephys
% 
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
fileList                 = {SpikeFileList; CaImagingShortDelayFastFileList; CaImagingShortDelaySlowFileList;...
                            CaImagingLongDelayFastFileList; CaImagingLongDelaySlowFileList; CaImagingShortDelaySlowVirusFileList};

nData                    = 1;
load([TempDatDir DataSetList(nData).name '.mat']);
ephysCellIndex{nData}        = [DataSetList(nData).cellinfo(:).depth]  > 110 &...
                           [DataSetList(nData).cellinfo(:).depth]  < 750;                        
for nData           = 2:length(fileList)
    load([TempDatDir DataSetList(nData).name '.mat']);
    ephysCellIndex{nData}  = [DataSetList(nData).cellinfo(:).AP_axis]  > 2100 &...
                      [DataSetList(nData).cellinfo(:).AP_axis]  < 2900 &...
                      [DataSetList(nData).cellinfo(:).ML_axis]  > 1100 &...
                      [DataSetList(nData).cellinfo(:).ML_axis]  < 1900;                        
end

save ([TempDatDir 'DataListShuffle.mat'], 'ephysCellIndex', 'fileList', '-append');
