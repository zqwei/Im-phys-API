%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
load('caData.mat')


if ~exist([PlotDir 'SingleUnitsTscore'],'dir')
    mkdir([PlotDir 'SingleUnitsTscore'])
end

cmap = cbrewer('qual', 'Set1', 9, 'cubic');
cmap = cmap([3, 5, 9], :);
groupColors = {cmap(3, :), cmap(2, :), cmap(1, :)};

load([TempDatDir DataSetList(1).name '.mat'])
unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(1).params);

numFold     = 30;

mono_       = [];
multi_      = [];

for frThres = [0 1 4 10] % spike count in this case
    load(['validMat_' num2str(frThres, '%02d')], 'validMat')
    sizeGroup = nan(numFold, 3);
    for nFold = 1:numFold
        nGroup = histcounts(unitGroup(validMat(:,nFold)), 0:3);
        sizeGroup(nFold, :) = nGroup(1:3)/sum(nGroup(1:3));
    end
    mono_   = [mono_, sizeGroup(:,2)];
    multi_  = [multi_, sizeGroup(:,3)];
end

mono_   = [mono_, caSizeGroup(:,2)];
multi_  = [multi_, caSizeGroup(:,3)];

figure
boxplot([mono_, multi_]);
ylim([0 0.7])
xlabel('Cell type')
ylabel('Frac. cell')
set(gca, 'TickDir', 'out')
setPrint(8, 6, 'Cell_type', 'pdf')
