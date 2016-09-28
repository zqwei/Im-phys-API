%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListC2SModel.mat']);

% cmap = cbrewer('qual', 'Set1', 3, 'cubic');
cmap = [0.8000    0.8000    0.8000;
       1.0000    0.6000         0;
       0    0.8000         0];

for nData      = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);    
    sizeGroup = histcounts(unitGroup, 0:3);
%     sizeGroup
%     sum(sizeGroup)
%     disp(sizeGroup(2)/sizeGroup(3))
    figure('Visible', 'off');
    groupNames      = {'Non.', 'Homo.', 'Dynamical'};
    pie(sizeGroup)
    colormap(cmap)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTime_' DataSetList(nData).name '_withOLRemoval'])
end
