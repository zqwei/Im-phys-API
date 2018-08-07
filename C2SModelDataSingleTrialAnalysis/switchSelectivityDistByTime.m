%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListC2SMCMCSingleTrialModel.mat']);

cmap = [0.8000    0.8000    0.8000;
       1.0000    0.6000         0;
       0    0.8000         0];

for nData      = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);    
    sizeGroup = histcounts(unitGroup, 0:3);
    figure('Visible', 'off');
    groupNames      = {'Non.', 'Homo.', 'Dynamical'};
    pie(sizeGroup)
    colormap(cmap)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTime_' DataSetList(nData).name '_withOLRemoval'])
end

load ([TempDatDir 'DataListShuffle.mat']);
nData = 10;
load([TempDatDir 'ModelSpikeMCMCSingleTrial_OOPSI_Ca_Fast_SShort_Delay.mat'])
unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
sizeGroup = histcounts(unitGroup, 0:3);
figure('Visible', 'off');
groupNames      = {'Non.', 'Homo.', 'Dynamical'};
pie(sizeGroup)
colormap(cmap)
set(gca, 'TickDir', 'out')
setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTime_ModelSpikeMCMCSingleTrial_OOPSI_Ca_Fast_SShort_Delay'])