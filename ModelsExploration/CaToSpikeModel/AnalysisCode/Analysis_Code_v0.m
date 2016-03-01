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

%  GCaMP5,6s,6f in vivo imaging/ephys data published in 
%  (Chen et. al. 2013 Nature; Akerboom, Chen 2012 J. Neurosci)

% 
% obtain the spike dataset from a list of files
% 
% version 1.0
%
% Comparison list
%
% Output:
% SpikeDataSet     --- yDim x 1 cells (yDims number of neurons) 
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Analysis_Code_v0 % generate the data of analysis

    addpath('../Func/')
    setDir;


    totCell     = repmat(struct('cellName',1, 'nRep', 1, 'expression', 'virus',...
                                'CaIndicator', 'GCaMP6f', 'fROI', 1, ...
                                'fNeuropil', 1, 'rawEphys', 1, 'filteredEphys',1, ...
                                'detectedSpikes', 1, 'CaTime', 1, 'ephysTime', 1),1000, 1);

    numCell     = 0;

    for nFolder = 1:length(dataFolders) %#ok<USENS>

        fileList     = dir([dataDir, dataFolders{nFolder}, '*.mat']);        
        for nCell    = 1:length(fileList)
           numCell                         = numCell + 1;
           fileName                        = fileList(nCell).name;
           [~, fileName, ~]                = fileparts(fileName);
           totCell(numCell).cellName       = fileName(1:end-4);
           totCell(numCell).nRep           = str2double(fileName(end-2:end));
           totCell(numCell).expression     = expression{nFolder}; %#ok<USENS>
           totCell(numCell).CaIndicator    = CaIndicator{nFolder}; %#ok<USENS>
           load([dataDir, dataFolders{nFolder}, fileList(nCell).name]);
           totCell(numCell).fROI           = obj.timeSeriesArrayHash.value{1}.valueMatrix;
           totCell(numCell).fNeuropil      = obj.timeSeriesArrayHash.value{2}.valueMatrix;
           totCell(numCell).rawEphys       = obj.timeSeriesArrayHash.value{3}.valueMatrix;
           totCell(numCell).filteredEphys  = obj.timeSeriesArrayHash.value{4}.valueMatrix;
           totCell(numCell).detectedSpikes = obj.timeSeriesArrayHash.value{5}.valueMatrix;
           totCell(numCell).CaTime         = obj.timeSeriesArrayHash.value{2}.time;
           totCell(numCell).ephysTime      = obj.timeSeriesArrayHash.value{5}.time;
        end        
    end

    totCell     = totCell(1:numCell); %#ok<NASGU>

    save([TempDir 'TotCell.mat'], 'totCell','-v7.3');