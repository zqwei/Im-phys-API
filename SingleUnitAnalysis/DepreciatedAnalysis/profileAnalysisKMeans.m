%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% profileAnalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function profileAnalysisKMeans
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
k             = 9;

if ~exist([PlotDir 'SingleUnitsActivityProfileKMeans'],'dir')
    mkdir([PlotDir 'SingleUnitsActivityProfileKMeans'])
end

for nData     = [1 3 4]
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotKMeans(nDataSet, k, DataSetList(nData).params);
    setPrint(8*ceil(sqrt(k)), 6*ceil(sqrt(k)), ...
        [PlotDir 'SingleUnitsActivityProfileKMeans/CentriodActivityProfileKMeans_K_' ...
        num2str(k, '%03d') '_' DataSetList(nData).name])
end

close all
