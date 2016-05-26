addpath('../Func');
setDir;

% if ~exist([PlotDir 'SingleUnitsRampingDown'],'dir')
%     mkdir([PlotDir 'SingleUnitsRampingDown'])
% end    
% 
% load ([TempDatDir 'DataListShuffle.mat']);
% nData = 1; % plot raster and psth
% load([TempDatDir DataSetList(nData).name '.mat'])
% params = DataSetList(nData).params;
% depth = [DataSetList(nData).cellinfo.depth];
% validDepth = depth<700 & depth>100;
% nDataSet   = nDataSet(validDepth);
% spikeDataSet = nDataSet;   
% 
% yesRampDown          = zeros(length(nDataSet), 1);
% noRampDown           = zeros(length(nDataSet), 1);
% timePoints(1)        = sum(DataSetList(nData).params.timeSeries<DataSetList(nData).params.polein);
% timePoints(2)        = sum(DataSetList(nData).params.timeSeries<0);
% contraIndex = false(length(nDataSet), 1);
% logPValueEpoch= getLogPValueSpikeEpoch(nDataSet, DataSetList(nData).params);
% unitGroup = getAnovaLogPSpikeEpoch(logPValueEpoch);
% 
% for nUnit            = 1:length(spikeDataSet)
%     meanPreSample    = [mean(spikeDataSet(nUnit).unit_yes_trial(:, 1:timePoints(1)), 2); mean(spikeDataSet(nUnit).unit_no_trial(:, 1:timePoints(1)), 2)];
%     meanYesSample    = mean(spikeDataSet(nUnit).unit_yes_trial(:, timePoints(1):timePoints(2)), 2);
%     meanNoSample     = mean(spikeDataSet(nUnit).unit_no_trial(:, timePoints(1):timePoints(2)), 2);
%     if ttest2(meanPreSample, meanYesSample, 'tail', 'right')
%         yesRampDown(nUnit) = 1;
%     elseif ttest2(meanPreSample, meanYesSample, 'tail', 'left')
%         yesRampDown(nUnit) = -1;
%     end
% 
%     if ttest2(meanPreSample, meanNoSample, 'tail', 'right')
%         noRampDown(nUnit)  = 1;
%     elseif ttest2(meanPreSample, meanNoSample, 'tail', 'left')
%         noRampDown(nUnit)  = -1;
%     end
%     
%     yesTrial = mean(nDataSet(nUnit).unit_yes_trial);
%     noTrial  = mean(nDataSet(nUnit).unit_no_trial);
%     contraIndex(nUnit)   = sum(noTrial(timePoints(2):end))<sum(yesTrial(timePoints(2):end));
% end
% 
% sum_yesRampDown = sum((yesRampDown==1 & noRampDown>=0) | (yesRampDown>=0 & noRampDown==1) & unitGroup>0);
% sum((yesRampDown==1 & noRampDown>=0) | (yesRampDown>=0 & noRampDown==1) & contraIndex & unitGroup>0)/sum_yesRampDown
% sum_noRampDown = sum((yesRampDown==-1 & noRampDown<=0) | (yesRampDown<=0 & noRampDown==-1)  & unitGroup>0);
% sum((yesRampDown==-1 & noRampDown<=0) | (yesRampDown<=0 & noRampDown==-1) & contraIndex & unitGroup>0)/sum_noRampDown
% 
% load ([TempDatDir 'DataListS2CModel.mat']);
% 
% nData = 2;
% load([TempDatDir DataSetList(nData).name '.mat'])
% nDataSet   = nDataSet(validDepth);
% logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
% unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);
% params      = DataSetList(nData).params;
% yesActMat   = nan(length(nDataSet), length(params.timeSeries));
% noActMat    = nan(length(nDataSet), length(params.timeSeries));
% timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
% 
% for nUnit   = 1:length(nDataSet)
%     yesTrial = mean(nDataSet(nUnit).unit_yes_trial);
%     noTrial  = mean(nDataSet(nUnit).unit_no_trial);
%     yesActMat(nUnit, :)  = yesTrial;
%     noActMat(nUnit, :)   = noTrial;
%     contraIndex(nUnit)   = sum(noTrial(timePoints(2):end))>sum(yesTrial(timePoints(2):end));
% end  
% 
% sum((yesRampDown==1 & noRampDown>=0) | (yesRampDown>=0 & noRampDown==1) & contraIndex & unitGroup>0)/sum_yesRampDown
% sum((yesRampDown==-1 & noRampDown<=0) | (yesRampDown<=0 & noRampDown==-1) & contraIndex & unitGroup>0)/sum_noRampDown


contraFrac = [0.4364, 0.5513];
selectiveFrac = [0.6910, 0.9289];
figure;
hold on
bar(1:2,[contraFrac.*selectiveFrac; contraFrac.*(1-selectiveFrac)]','stack')
gridxy([], [0.5], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
xlim([0.5 2.5])
ylim([0 0.65])
xlabel('Dataset index')
ylabel('Frac. Contra')
box off;
setPrint(4, 3, [PlotDir 'S2CModel/ContraFractionRampDownUpLongDecay'])


contraFrac = [0.4364, 0.5513];
selectiveFrac = [0.8134, 0.9289];
figure;
hold on
bar(1:2,[contraFrac.*selectiveFrac; contraFrac.*(1-selectiveFrac)]','stack')
gridxy([], [0.5], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
xlim([0.5 2.5])
ylim([0 0.65])
xlabel('Dataset index')
ylabel('Frac. Contra')
box off;
setPrint(4, 3, [PlotDir 'S2CModel/ContraFractionRampDownUpShortDecay'])
