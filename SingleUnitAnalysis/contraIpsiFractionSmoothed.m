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

load([TempDatDir DataSetList(nData).name '_old.mat'])
depth = [nDataSet.depth_in_um];
validDepth = depth<800 & depth>100;

% selective defined as neuron shows selectivity for any of sample, delay,
% response epoch.
logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
unitGroup   = plotTtestLogPSpikeEpoch (logPValueEpoch)'; 
params      = DataSetList(nData).params;
contraIndex = false(length(nDataSet), 1);
cellType    = [nDataSet.cell_type];
yesActMat   = nan(length(nDataSet), length(smoothedparams.timeSeries));
noActMat    = nan(length(nDataSet), length(smoothedparams.timeSeries));
timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);

% contra vs ipsi as neuron firing difference summed across all thress
% epochs
for nUnit   = 1:length(nDataSet)
    yesTrial = mean(nDataSet(nUnit).unit_yes_trial);
    noTrial  = mean(nDataSet(nUnit).unit_no_trial);
    yesActMat(nUnit, :)  = mean(smoothDataSet(nUnit).unit_yes_trial)/1000;
    noActMat(nUnit, :)   = mean(smoothDataSet(nUnit).unit_no_trial)/1000;
    contraIndex(nUnit)   = sum(noTrial(timePoints(2):end))<sum(yesTrial(timePoints(2):end));
end

otherIndex = unitGroup~=0 & cellType == 1 & validDepth;
otherIndex = otherIndex';
    
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
set(gca, 'TickDir', 'out')

actPre = mean([yesActMat(contraIndex& otherIndex,1:100), noActMat(contraIndex & otherIndex,1:100)], 2);
actYes = mean(yesActMat(contraIndex& otherIndex,435:550), 2);
actNo = mean(noActMat(contraIndex& otherIndex,435:550), 2);

[~, p] = ttest2(actPre, actYes, 'tail', 'left');
[~, p] = ttest(actPre, actNo, 'tail', 'right');

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
set(gca, 'TickDir', 'out')


actPre = mean([yesActMat(~contraIndex& otherIndex,1:100), noActMat(~contraIndex & otherIndex,1:100)], 2);
actYes = mean(yesActMat(~contraIndex& otherIndex,435:550), 2);
actNo = mean(noActMat(~contraIndex& otherIndex,435:550), 2);

[~, p] = ttest2(actPre, actYes, 'tail', 'right');
[~, p] = ttest(actPre, actNo, 'tail', 'left');

setPrint(8*2, 6, [PlotDir 'SingleUnitsContraIpsi\' DataSetList(nData).name])
