%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% no-selective cell      : no selectivity -- 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load([TempDatDir 'DataListEphys.mat'], 'DataSetList');
smoothedparams = DataSetList(5).params;
load ([TempDatDir 'DataListShuffle.mat']);


if ~exist([PlotDir 'SingleUnitsContraIpsi'],'dir')
    mkdir([PlotDir 'SingleUnitsContraIpsi'])
end

cmap = cbrewer('qual', 'Set1', 3, 'cubic');

nData = 1;
load([TempDatDir 'Shuffle_Spikes_boxcar_200ms.mat'])

smoothDataSet = nDataSet;

load([TempDatDir DataSetList(nData).name '.mat'])
depth = [DataSetList(nData).cellinfo.depth];
validDepth = depth<700 & depth>100;

logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);
params      = DataSetList(nData).params;
contraIndex = false(length(nDataSet), 1);
cellType    = [DataSetList(nData).cellinfo.cellType]';
yesActMat   = nan(length(nDataSet), length(smoothedparams.timeSeries));
noActMat    = nan(length(nDataSet), length(smoothedparams.timeSeries));
timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);

for nUnit   = 1:length(nDataSet)
    yesTrial = mean(nDataSet(nUnit).unit_yes_trial);
    noTrial  = mean(nDataSet(nUnit).unit_no_trial);
    yesActMat(nUnit, :)  = mean(smoothDataSet(nUnit).unit_yes_trial)/1000;
    noActMat(nUnit, :)   = mean(smoothDataSet(nUnit).unit_no_trial)/1000;
    contraIndex(nUnit)   = sum(noTrial(timePoints(2):end))<sum(yesTrial(timePoints(2):end));
end

cellType  = cellType_all(usedUnits==1,1);
otherIndex = unitGroup~=0 & cellType == 1; % & validDepth';
    
figure;
subplot(1, 2, 1)
hold on
shadedErrorBar(smoothedparams.timeSeries, mean(noActMat(contraIndex & otherIndex,:)), sem(noActMat(contraIndex& otherIndex,:)),'-r')
shadedErrorBar(smoothedparams.timeSeries, mean(yesActMat(contraIndex& otherIndex,:)), sem(yesActMat(contraIndex& otherIndex,:)),'-b')
ylim([2, 10])
title('contra neuron only')
gridxy ([smoothedparams.polein, smoothedparams.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
xlim([smoothedparams.timeSeries(2) smoothedparams.timeSeries(end-30)])
ylim([2, 10])
xlabel('Time (ms)')
box off
ylabel('Mean activity')

% baseline1 = mean(mean(noActMat(contraIndex, 1:8)));
% baseline2 = mean(mean(yesActMat(contraIndex, 1:8)));
% baseline  = (baseline1 + baseline2)/2;
% 
% disp(max(mean(noActMat(contraIndex,:))) - baseline)
% disp(min(mean(yesActMat(contraIndex,:))) - baseline)

subplot(1, 2, 2)
hold on
shadedErrorBar(smoothedparams.timeSeries, mean(noActMat(~contraIndex& otherIndex,:)), sem(noActMat(~contraIndex& otherIndex,:)),'-r')
shadedErrorBar(smoothedparams.timeSeries, mean(yesActMat(~contraIndex& otherIndex,:)), sem(yesActMat(~contraIndex& otherIndex,:)),'-b')
ylim([2, 10])
title('ipsi neuron only')
gridxy ([smoothedparams.polein, smoothedparams.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
xlim([smoothedparams.timeSeries(2) smoothedparams.timeSeries(end-30)])
ylim([2, 10])
xlabel('Time (ms)')
box off
ylabel('Mean activity')

% baseline1 = mean(mean(noActMat(~contraIndex, 1:8)));
% baseline2 = mean(mean(yesActMat(~contraIndex, 1:8)));
% baseline  = (baseline1 + baseline2)/2;
% 
% disp(min(mean(noActMat(~contraIndex,:))) - baseline)
% disp(max(mean(yesActMat(~contraIndex,:))) - baseline)

setPrint(8*2, 6, [PlotDir 'SingleUnitsContraIpsi\' DataSetList(nData).name])
