%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selectiviy index of Nuo's version
% Nomarilized contra - ispi firing rate
% 
% For this analysis, firing rate data for spiking neurons are recomputed
% using Nuo's filter
% 
% Those for Ca++ remain the same as the other analysis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);

%%%% For spiking data

if ~exist([PlotDir 'ModeledComparingSingleUnitsNormalizedDiff'],'dir')
    mkdir([PlotDir 'ModeledComparingSingleUnitsNormalizedDiff'])
end

params.frameRate       =  1000;
params.binsize         =  1/params.frameRate; % 1 ms bin
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round(-3.5 * params.frameRate); % these parameters from Nuo
maxTimeToAnalysis      =  round(2.0 * params.frameRate); % these parameters from Nuo
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  20;
params.expression      = 'None';
nDataSet               = getSpikeData(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);

cellType               = false(length(nDataSet), 1);

for nUnit  = 1:length(nDataSet)
    if strcmp(nDataSet(nUnit).cell_type, 'putative_interneuron')
        cellType(nUnit)   = false;
    elseif strcmp(nDataSet(nUnit).cell_type, 'putative_pyramidal')
        cellType(nUnit)   = true;
    end   
end

pryCellIndex    = cellType;
nDataSet        = nDataSet(pryCellIndex);
PSTHStartTime   = -3.5;
PSTHEndTime     = 2.2;
wholeTrialStartingTime = -3.1; % pole in = -2.6 % pole out = -1.4
wholeTrialEndTime      = 1.3; % these two values are from Nuo

psthTime             = PSTHStartTime:.001:PSTHEndTime;
psthTime             = psthTime(201:end-200);
periodIndex(1)       = sum(psthTime<wholeTrialStartingTime);
periodIndex(2)       = sum(psthTime<params.polein);
periodIndex(3)       = sum(psthTime<params.poleout);
periodIndex(4)       = sum(psthTime<0);
periodIndex(5)       = sum(psthTime<wholeTrialEndTime);
sigSelective         = cell2mat(arrayfun(@(x) getSelective(x, periodIndex),...
                                        nDataSet, 'UniformOutput', false));

cellType             = ones(length(nDataSet),1) * 4;
cellType( (sigSelective(:,1)|sigSelective(:,2))  & (~sigSelective(:,3)) ) = 1;
cellType( (sigSelective(:,1)|sigSelective(:,2))  &  sigSelective(:,3))    = 2;
cellType( ~(sigSelective(:,1)|sigSelective(:,2))  & sigSelective(:,3))    = 3;
spkCellType          = cellType;


for nData = 2:length(DataSetList)    
    load([TempDatDir DataSetList(nData).name '.mat'])
    nDataSet             = nDataSet(pryCellIndex);
    psthTime             = DataSetList(nData).params.timeSeries;
    periodIndex(1)       = sum(psthTime<wholeTrialStartingTime);
    periodIndex(2)       = sum(psthTime<params.polein);
    periodIndex(3)       = sum(psthTime<params.poleout);
    periodIndex(4)       = sum(psthTime<0);
    periodIndex(5)       = sum(psthTime<wholeTrialEndTime);
    sigSelective         = cell2mat(arrayfun(@(x) getSelective(x, periodIndex),...
                                        nDataSet, 'UniformOutput', false));
    cellType             = ones(length(nDataSet),1) * 4;
    cellType( (sigSelective(:,1)|sigSelective(:,2))  & (~sigSelective(:,3)) ) = 1;
    cellType( (sigSelective(:,1)|sigSelective(:,2))  &  sigSelective(:,3))    = 2;
    cellType( ~(sigSelective(:,1)|sigSelective(:,2))  & sigSelective(:,3))    = 3;
    
    groupCounts = zeros(4, 4);
    
    for nGroup  = 1:4
        for mGroup = 1:4
            groupCounts(nGroup, mGroup) = sum(cellType == mGroup & spkCellType == nGroup);
        end
    end
    
    groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
    
    figure;
    barh(1:4, groupPerCounts, 'stack', 'edgecolor', 'none');
    caxis([1 4])
    xlim([0 1])
    box off
    xlabel('% cell type')
    ylabel('cell type')
    ylim([0.5 4.5])
    set(gca, 'yTick', 1:8, 'yTickLabel', {'Pre', 'Peri+Pre', 'Peri', 'Non'})
    setPrint(8, 6, [PlotDir 'ModeledComparingSingleUnitsNormalizedDiff/SingleUnitsNormalizedDiff_' DataSetList(nData).name], 'pdf')
end

close all