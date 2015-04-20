% 
% Comparison based on single unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0.2
% 
% 
% 
% Gather all information of each analyzed cell
% print out the information like depth range, number of trials, number of
% cells
% 
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);


fileList            = {SpikeFileList; CaImagingShortDelayFastFileList; CaImagingShortDelaySlowFileList;...
                       CaImagingLongDelayFastFileList; CaImagingLongDelaySlowFileList; CaImagingShortDelaySlowVirusFileList};
                   
nData               = 1;
load([TempDatDir DataSetList(nData).name '.mat']);
pyrCellIndex        = [DataSetList(nData).cellinfo(:).cellType]   == 1;  
inhCellIndex        = [DataSetList(nData).cellinfo(:).cellType]   == 0; 
disp('num of pyrCell:')
disp(sum(pyrCellIndex))
disp('num of inhCell:')
disp(sum(inhCellIndex))
disp('pyrCellDepthRange:')
disp(num2str(min([DataSetList(nData).cellinfo(pyrCellIndex).depth]), '%.3f')) 
disp(num2str(max([DataSetList(nData).cellinfo(pyrCellIndex).depth]), '%.3f')) 
disp('inhCellDepthRange:')
disp(num2str(min([DataSetList(nData).cellinfo(inhCellIndex).depth]), '%.3f')) 
disp(num2str(max([DataSetList(nData).cellinfo(inhCellIndex).depth]), '%.3f')) 
numOfYesFile        = cell2mat(arrayfun(@(x) length(x.unit_yes_trial_index), nDataSet, 'UniformOutput', false));
numOfNoFile         = cell2mat(arrayfun(@(x) length(x.unit_no_trial_index), nDataSet, 'UniformOutput', false));
disp('pyrCellYesFile:')
disp(mean(numOfYesFile(pyrCellIndex)))
disp('pyrCellNoFile:')
disp(mean(numOfNoFile(pyrCellIndex)))
disp('inhCellYesFile:')
disp(mean(numOfYesFile(inhCellIndex)))
disp('inhCellNoFile:')
disp(mean(numOfNoFile(inhCellIndex)))

for nData           = 2:length(fileList)
    disp(DataSetList(nData).name)
    load([TempDatDir DataSetList(nData).name '.mat']);
    ephysCellIndex  = [DataSetList(nData).cellinfo(:).AP_axis]  > 2100 &...
                      [DataSetList(nData).cellinfo(:).AP_axis]  < 2900 &...
                      [DataSetList(nData).cellinfo(:).ML_axis]  > 1100 &...
                      [DataSetList(nData).cellinfo(:).ML_axis]  < 1900;  
    otherCellIndex  = ~ephysCellIndex; 
    disp('num of ephysCell:')
    disp(sum(ephysCellIndex))
    disp('num of otherCell:')
    disp(sum(otherCellIndex))
    disp('ephysCellDepthRange:')
    disp(num2str(min([DataSetList(nData).cellinfo(ephysCellIndex).depth]), '%.3f')) 
    disp(num2str(max([DataSetList(nData).cellinfo(ephysCellIndex).depth]), '%.3f')) 
    disp('otherCellDepthRange:')
    disp(num2str(min([DataSetList(nData).cellinfo(otherCellIndex).depth]), '%.3f')) 
    disp(num2str(max([DataSetList(nData).cellinfo(otherCellIndex).depth]), '%.3f')) 
    disp('otherCellAPRange:')
    disp(num2str(min([DataSetList(nData).cellinfo(otherCellIndex).AP_axis]/1000), '%.3f')) 
    disp(num2str(max([DataSetList(nData).cellinfo(otherCellIndex).AP_axis]/1000), '%.3f')) 
    disp('otherCellMLRange:')
    disp(num2str(min([DataSetList(nData).cellinfo(otherCellIndex).ML_axis]/1000), '%.3f')) 
    disp(num2str(max([DataSetList(nData).cellinfo(otherCellIndex).ML_axis]/1000), '%.3f')) 
    numOfYesFile    = cell2mat(arrayfun(@(x) length(x.unit_yes_trial_index), nDataSet, 'UniformOutput', false));
    numOfNoFile     = cell2mat(arrayfun(@(x) length(x.unit_no_trial_index), nDataSet, 'UniformOutput', false));
    disp('ephysCellYesFile:')
    disp(mean(numOfYesFile(ephysCellIndex)))
    disp('ephysCellNoFile:')
    disp(mean(numOfNoFile(ephysCellIndex)))
    disp('otherCellYesFile:')
    disp(mean(numOfYesFile(otherCellIndex)))
    disp('otherCellNoFile:')
    disp(mean(numOfNoFile(otherCellIndex)))
end
