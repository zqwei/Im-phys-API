addpath('../Func');
setDir;

if ~exist([PlotDir 'SingleUnitsRampingDown'],'dir')
    mkdir([PlotDir 'SingleUnitsRampingDown'])
end    

load([TempDatDir 'DataListShuffle.mat']);
nData = 1;
load([TempDatDir DataSetList(nData).name '.mat'])
depth       = [nDataSet.depth_in_um];
spikeDataSet= nDataSet(depth<471);  

yesRampDown          = zeros(length(spikeDataSet), 1);
noRampDown           = zeros(length(spikeDataSet), 1);
timePoints(1)        = sum(DataSetList(nData).params.timeSeries<DataSetList(nData).params.polein);
timePoints(2)        = sum(DataSetList(nData).params.timeSeries<0);
contraIndex          = false(length(spikeDataSet), 1);
unitGroup            = getLogPValueTscoreSpikeTime(spikeDataSet, DataSetList(nData).params); 

for nUnit            = 1:length(spikeDataSet)
    meanPreSample    = [mean(spikeDataSet(nUnit).unit_yes_trial(:, 1:timePoints(1)), 2); mean(spikeDataSet(nUnit).unit_no_trial(:, 1:timePoints(1)), 2)];
    meanYesSample    = mean(spikeDataSet(nUnit).unit_yes_trial(:, timePoints(1):timePoints(2)), 2);
    meanNoSample     = mean(spikeDataSet(nUnit).unit_no_trial(:, timePoints(1):timePoints(2)), 2);
    if ttest2(meanPreSample, meanYesSample, 'tail', 'right')
        yesRampDown(nUnit) = 1;
    elseif ttest2(meanPreSample, meanYesSample, 'tail', 'left')
        yesRampDown(nUnit) = -1;
    end    
%     yesRampDown(nUnit) = mean(meanYesSample)<mean(meanPreSample);
    if ttest2(meanPreSample, meanNoSample, 'tail', 'right')
        noRampDown(nUnit)  = 1;
    elseif ttest2(meanPreSample, meanNoSample, 'tail', 'left')
        noRampDown(nUnit)  = -1;
    end    
%     noRampDown(nUnit) = mean(meanNoSample)<mean(meanPreSample);
    yesTrial = mean(nDataSet(nUnit).unit_yes_trial);
    noTrial  = mean(nDataSet(nUnit).unit_no_trial);
    contraIndex(nUnit)= sum(noTrial(timePoints(1):end))<sum(yesTrial(timePoints(1):end));
end


RampDown = (yesRampDown==1 & noRampDown>=0) | (yesRampDown>=0 & noRampDown==1) & unitGroup>0;
RampUp   = (yesRampDown==-1 & noRampDown<=0) | (yesRampDown<=0 & noRampDown==-1)  & unitGroup>0;

spkContraIndex = contraIndex;
% sum(RampDown)
% sum(RampUp)
% sum(RampDown & contraIndex)
% sum(RampUp & contraIndex)


% [~, ~, anmIndex] = unique(cell2mat({DataSetList(nData).cellinfo.anmName}'), 'rows');
% numGroup    = 3;
% groupCounts = zeros(anmIndex(end), numGroup);
% for nAnm    = 1:anmIndex(end)
%     groupCounts(nAnm, 1) = sum(RampDown & anmIndex == nAnm & contraIndex);
%     groupCounts(nAnm, 2) = sum(RampUp & anmIndex == nAnm & contraIndex);
%     groupCounts(nAnm, 3) = sum(unitGroup>0 & anmIndex == nAnm & contraIndex) - sum(groupCounts(nAnm, [1 2]));
%     groupCounts(nAnm, 3) = max(groupCounts(nAnm, 3), 0);
%     groupCounts(nAnm, :) = groupCounts(nAnm, :)/sum(groupCounts(nAnm, :));
% end
% 
% groupContra = groupCounts;
% 
% for nAnm    = 1:anmIndex(end)
%     groupCounts(nAnm, 1) = sum(RampDown & anmIndex == nAnm & ~contraIndex);
%     groupCounts(nAnm, 2) = sum(RampUp & anmIndex == nAnm & ~contraIndex);
%     groupCounts(nAnm, 3) = sum(unitGroup>0 & anmIndex == nAnm & ~contraIndex) - sum(groupCounts(nAnm, [1 2]));
%     groupCounts(nAnm, 3) = max(groupCounts(nAnm, 3), 0);
%     groupCounts(nAnm, :) = groupCounts(nAnm, :)/sum(groupCounts(nAnm, :));
% end
% 
% groupIpsi = groupCounts;



load ([TempDatDir 'DataListS2C6fModel.mat']);

nData = 2;
load([TempDatDir DataSetList(nData).name '.mat'])
unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params); 

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
%     contraIndex(nUnit)   = sum(noTrial(timePoints(1):timePoints(2)))>sum(yesTrial(timePoints(1):timePoints(2)));
% end  

sum(RampDown & spkContraIndex & unitGroup)
sum(RampDown & ~spkContraIndex  & unitGroup)
sum(RampUp & spkContraIndex & unitGroup)
sum(RampUp & ~spkContraIndex  & unitGroup)

% first Down, second Up
contraFrac = [6/19, 20/34];
ipsiFrac = [9/23, 17/38];
figure;
hold on
errorbar(1:2,contraFrac, sqrt(contraFrac.*(1-contraFrac)/42),'-or','markerfacecolor','r')
errorbar(1:2,ipsiFrac, sqrt(ipsiFrac.*(1-ipsiFrac)/72),'-ob','markerfacecolor','b')
hold off;
xlim([0.5 2.5])
ylim([0.2 0.9])
xlabel('Dataset index')
ylabel('Frac. Contra')
box off;
set(gca, 'TickDir', 'out')
setPrint(6, 6, [PlotDir 'SingleUnitsRampingDown/ContraFractionRampDownUp6f'])
% 
% 
% contraFrac = [35/76, 110/140];
% ipsiFrac = [53/110, 81/100];
% figure;
% hold on
% errorbar(1:2,contraFrac, sqrt(contraFrac.*(1-contraFrac)/38),'-or','markerfacecolor','r')
% errorbar(1:2,ipsiFrac, sqrt(ipsiFrac.*(1-ipsiFrac)/38),'-ob','markerfacecolor','b')
% hold off;
% xlim([0.5 2.5])
% ylim([0.3 0.9])
% xlabel('Dataset index')
% ylabel('Frac. Contra')
% set(gca, 'TickDir', 'out')
% box off;
% setPrint(8, 6, [PlotDir 'SingleUnitsRampingDown/ContraFractionRampDownUpShortDecay'])
% 
% % down, up, other
% contraFrac = [76 140 90]/306;
% ipsiFrac = [110 100 125]/335;
% contraStd = sem(groupContra);
% ipsiStd = sem(groupIpsi);
% 
% figure;
% hold on
% plot(1:3, groupContra, '-o', 'color', [0.5 0.5 0.5])
% bar(1:3,contraFrac, 'edgecolor', 'k', 'facecolor', 'none')
% errorbar(1:3, contraFrac, contraStd, '.k')
% hold off;
% xlim([0.5 3.5])
% ylim([0 1])
% xlabel('Dataset index')
% ylabel('Frac. Contra')
% set(gca, 'TickDir', 'out')
% box off;
% setPrint(8, 6, [PlotDir 'SingleUnitsRampingDown/ContraFractionNumNeuron'])
% 
% figure;
% hold on
% plot(1:3, groupIpsi, '-o', 'color', [0.5 0.5 0.5])
% bar(1:3,ipsiFrac, 'edgecolor', 'k', 'facecolor', 'none')
% errorbar(1:3, ipsiFrac, contraStd, '.k')
% hold off;
% xlim([0.5 3.5])
% ylim([0 1])
% xlabel('Dataset index')
% ylabel('Frac. Contra')
% set(gca, 'TickDir', 'out')
% box off;
% setPrint(8, 6, [PlotDir 'SingleUnitsRampingDown/IpsiFractionNumNeuron'])
