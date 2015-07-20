% 
% Comparison based on single unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0.1
%


% Gather all information of each analyzed cell

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

% -------------------------------------------------------------------------
% Spiking Dataset
% -------------------------------------------------------------------------

fileList            = {SpikeFileList; CaImagingShortDelayFastFileList; CaImagingShortDelaySlowFileList;...
                       CaImagingLongDelayFastFileList; CaImagingLongDelaySlowFileList; CaImagingShortDelaySlowVirusFileList};


for nData           = [1 3 6]%1:length(fileList)
    load([TempDatDir DataSetList(nData).name '.mat']);
    DataSetList(nData).cellinfo  = repmat(struct('fileName',1, 'nUnit', 1, ...
                                'AP_axis',1, 'ML_axis', 1, 'depth', 1,...
                                'expression','none', 'cellType',-1),length(nDataSet), 1);     %#ok<SAGROW>
    for nUnit  = 1:length(nDataSet)
        nFileList                                     = fileList{nData};
        DataSetList(nData).cellinfo(nUnit).fileName   = nFileList(nDataSet(nUnit).sessionIndex).name;
        DataSetList(nData).cellinfo(nUnit).nUnit      = nDataSet(nUnit).nUnit;
        DataSetList(nData).cellinfo(nUnit).AP_axis    = nDataSet(nUnit).AP_in_um;
        DataSetList(nData).cellinfo(nUnit).ML_axis    = nDataSet(nUnit).ML_in_um;
        DataSetList(nData).cellinfo(nUnit).depth      = nDataSet(nUnit).depth_in_um;
        DataSetList(nData).cellinfo(nUnit).expression = DataSetList(nData).params.expression;
        if strcmp(nDataSet(nUnit).cell_type, 'putative_interneuron')
            DataSetList(nData).cellinfo(nUnit).cellType   = 0;
        elseif strcmp(nDataSet(nUnit).cell_type, 'putative_pyramidal')
            DataSetList(nData).cellinfo(nUnit).cellType   = 1;
        end        
    end
end

save([TempDatDir 'DataListShuffle.mat'], 'DataSetList');