%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate data from simultaneous recording of V1 cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Collecting raw data from TCW's file
% Includes:
% 1. Cell Name
% 2. n Rep
% 3. Expression
% 4. Ca++ indicator
% 5. ROI
% 6. Neuropil
% 7. Raw Ephys
% 8. Filtered Ephys % high pass
% 9. Detected Spikes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GCaMP5,6s,6f in vivo imaging/ephys data published in 
%  (Chen et. al. 2013 Nature; Akerboom, Chen 2012 J. Neurosci)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% obtain the spike dataset from a list of files
% 
% version 1.0
%
% Comparison list
%
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% informations used for old version of code
% totCell     = repmat(struct('cellName',1, 'nRep', 1, 'expression', 'virus',...
%                             'CaIndicator', 'GCaMP6f', 'fROI', 1, ...
%                             'fNeuropil', 1, 'rawEphys', 1, 'filteredEphys',1, ...
%                             'detectedSpikes', 1, 'CaTime', 1, 'ephysTime', 1),1000, 1);
% 
% numCell     = 0;
% 
% for nFolder = 1:length(dataFolders)
% 
%     fileList     = dir([dataDir, dataFolders{nFolder}, '*.mat']);        
%     for nCell    = 1:length(fileList)
%        numCell                         = numCell + 1;
%        fileName                        = fileList(nCell).name;
%        [~, fileName, ~]                = fileparts(fileName);
%        totCell(numCell).cellName       = fileName(1:end-4);
%        totCell(numCell).nRep           = str2double(fileName(end-2:end));
%        totCell(numCell).expression     = expression{nFolder};
%        totCell(numCell).CaIndicator    = CaIndicator{nFolder};
%        load([dataDir, dataFolders{nFolder}, fileList(nCell).name]);
%        totCell(numCell).fROI           = obj.timeSeriesArrayHash.value{1}.valueMatrix;
%        totCell(numCell).fNeuropil      = obj.timeSeriesArrayHash.value{2}.valueMatrix;
%        totCell(numCell).rawEphys       = obj.timeSeriesArrayHash.value{3}.valueMatrix;
%        totCell(numCell).filteredEphys  = obj.timeSeriesArrayHash.value{4}.valueMatrix;
%        totCell(numCell).detectedSpikes = obj.timeSeriesArrayHash.value{5}.valueMatrix;
%        totCell(numCell).CaTime         = obj.timeSeriesArrayHash.value{2}.time;
%        totCell(numCell).ephysTime      = obj.timeSeriesArrayHash.value{5}.time;
%     end        
% end
% 
% totCell     = totCell(1:numCell);

% save([TempDir 'TotCell.mat'], 'totCell','-v7.3');

%% informations --
%                 DF/F for each cell
%                 spk time for each cell


addpath('../Func');
setDirV1Cells;


totCell     = repmat(struct('cellName',1, 'nRep', 1, 'expression', 'virus',...
                            'CaIndicator', 'GCaMP6f', 'spk', 1, 'dff', 1),1000, 1);

numCell     = 0;

for nFolder = 1:length(dataFolders)

fileList     = dir([dataDir, dataFolders{nFolder}, '*.mat']);        
    for nCell    = 1:length(fileList)        
        fileName             = fileList(nCell).name;
        [~, fileName, ~]     = fileparts(fileName);
        load([dataDir, dataFolders{nFolder}, fileList(nCell).name]);
        CaTime               = obj.timeSeriesArrayHash.value{2}.time;
        fROI                 = obj.timeSeriesArrayHash.value{1}.valueMatrix;
        fNeuropil            = obj.timeSeriesArrayHash.value{2}.valueMatrix;
        rawEphys             = obj.timeSeriesArrayHash.value{3}.valueMatrix;
        filteredEphys        = obj.timeSeriesArrayHash.value{4}.valueMatrix;
        detectedSpikes       = obj.timeSeriesArrayHash.value{5}.valueMatrix;
        ephysTime            = obj.timeSeriesArrayHash.value{5}.time;

        if length(CaTime) == 14400
            numCell                      = numCell + 1;
            totCell(numCell).cellName    = fileName(1:end-4);
            totCell(numCell).nRep        = str2double(fileName(end-2:end));
            totCell(numCell).expression  = expression{nFolder};
            totCell(numCell).CaIndicator = CaIndicator{nFolder};
            para.t_frame                 = CaTime;
            para.t_ephys                 = ephysTime;
            para.fneuropil               = fNeuropil;
            para.fmean                   = fROI;
            para.filt                    = filteredEphys;
            para.peak                    = detectedSpikes;
            dff                          = get_baseline_corr_dff(para);
            spk                          = para.t_ephys(para.peak);            
            totCell(numCell).CaTime      = CaTime;
            totCell(numCell).spk         = spk;
            totCell(numCell).dff         = dff;
        end
    end        
end

totCell     = totCell(1:numCell);

save([TempDatDir 'DataListCells.mat'], 'totCell');