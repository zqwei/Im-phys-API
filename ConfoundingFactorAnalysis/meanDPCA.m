%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffleConfounding.mat']);

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;


numComps             = 15;
if ~exist([PlotDir 'ConfoundingFactordPCA'],'dir')
    mkdir([PlotDir 'ConfoundingFactordPCA'])
end

DataSetToAnalysis = [1 5 6];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All area vs area of spiking recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nData          = 5;
load([TempDatDir DataSetList(nData).name '.mat']);
matFileName1   = [TempDatDir 'optimalLambdas_' DataSetList(nData).name '_ConfoundingFactorAllData.mat'];
matFileName2   = [TempDatDir 'ConfoundingFactorddPCA_' DataSetList(nData).name '_AllData.mat'];
figFileName1   = [PlotDir 'ConfoundingFactordPCA/CollectedUnitsdPCAPMLAll_' DataSetList(nData).name ];
figFileName2   = [PlotDir 'ConfoundingFactordPCA/CollectedUnitsPCAAPMLAll_' DataSetList(nData).name ];
plotdPCAEV(nDataSet, numTrials, numComps, combinedParams, matFileName1, matFileName2, figFileName1, figFileName2)

APLoc          = [DataSetList(nData).cellinfo.AP_axis];
MLLoc          = [DataSetList(nData).cellinfo.ML_axis];
nDataSet       = nDataSet(APLoc>2400 & APLoc<2600 & MLLoc>1100 & MLLoc<1900);
matFileName1   = [TempDatDir 'optimalLambdas_' DataSetList(nData).name '_ConfoundingFactorSubData.mat'];
matFileName2   = [TempDatDir 'ConfoundingFactorddPCA_' DataSetList(nData).name '_SubData.mat'];
figFileName1   = [PlotDir 'ConfoundingFactordPCA/CollectedUnitsdPCAPMLSub_' DataSetList(nData).name ];
figFileName2   = [PlotDir 'ConfoundingFactordPCA/CollectedUnitsPCAAPMLSub_' DataSetList(nData).name ];
plotdPCAEV(nDataSet, numTrials, numComps, combinedParams, matFileName1, matFileName2, figFileName1, figFileName2)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison across animals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for nData              = DataSetToAnalysis
    load([TempDatDir DataSetList(nData).name '.mat']);
    [~, ~, anmIndex] = unique(cell2mat({DataSetList(nData).cellinfo.anmName}'), 'rows');
    
    for nAnm    = 1:anmIndex(end)
        if sum(anmIndex == nAnm) > 50
            tDataSet = nDataSet(anmIndex == nAnm);
            matFileName1   = [TempDatDir 'optimalLambdas_' DataSetList(nData).name '_ANM_' num2str(nAnm) '.mat'];
            matFileName2   = [TempDatDir 'ConfoundingFactorddPCA_' DataSetList(nData).name '_ANM_' num2str(nAnm) '.mat'];
            figFileName1   = [PlotDir 'ConfoundingFactordPCA/CollectedUnitsdPCAANM_' DataSetList(nData).name  '_ANM_' num2str(nAnm) ];
            figFileName2   = [PlotDir 'ConfoundingFactordPCA/CollectedUnitsPCAANM_' DataSetList(nData).name  '_ANM_' num2str(nAnm) ];
            plotdPCAEV(tDataSet, numTrials, numComps, combinedParams, matFileName1, matFileName2, figFileName1, figFileName2)
        end
    end
    
    
end

close all