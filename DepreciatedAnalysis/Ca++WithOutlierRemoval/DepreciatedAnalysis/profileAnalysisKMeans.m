%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% profileAnalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function profileAnalysisKMeans
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
k             = 9;

for nData     = [3 4]
    load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat']) %
    plotKMeans(nDataSet, k, DataSetList(nData).params);
    setPrint(8*ceil(sqrt(k)), 6*ceil(sqrt(k)), [PlotDir 'SingleUnitsActivityProfileKMeans/CentriodActivityProfileKMeans_K_' num2str(k, '%03d') '_' DataSetList(nData).name  '_withOLRemoval'])
end

close all