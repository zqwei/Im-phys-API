%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);


if ~exist([PlotDir 'SingleUnitsTscore'],'dir')
    mkdir([PlotDir 'SingleUnitsTscore'])
end


nData   = 1;
load([TempDatDir DataSetList(nData).name '.mat'])
unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params); 
fr = nan(length(nDataSet), 1);
for nUnit = 1:length(nDataSet)
    meanYes = mean(mean(nDataSet(nUnit).unit_yes_trial));
    meanNo  = mean(mean(nDataSet(nUnit).unit_no_trial));
    fr(nUnit) = (meanYes + meanNo)/2;
end

% boxplot(fr, unitGroup)
% box off
% ylabel('Mean firing rate')
% xlabel('Cell group')
% set(gca, 'TickDir', 'out')
% setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTimeFiringRate_' DataSetList(nData).name])
% 
