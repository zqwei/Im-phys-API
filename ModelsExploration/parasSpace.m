%
% parameter space for neuronal trace from S2C
%
%

addpath('../Func');
setDirV1Cells;
load([TempDatDir 'DataListCells.mat'], 'totCell');

if ~exist([PlotDir 'ModelCellFits'],'dir')
    mkdir([PlotDir 'ModelCellFits'])
end

load([TempDatDir 'ParamsFitCells.mat'], 'paras');

nKeys    = {'tau_r', 'tau_d', 'n', 'K', 'Fm'};

paraMat  = nan(length(paras), length(nKeys));

for nKey = 1:length(nKeys)
    paraMat(:, nKey) = [paras.(nKeys{nKey})];
end

paraMat(paraMat(:,1)>0.2, 1) = nan;
paraMat(paraMat(:,5)<0 | paraMat(:,5)>500, 5) = nan;
paraMat(paraMat(:,4)>200, 4) = nan;
paraMat(paraMat(:,3)>4, 3) = nan;

group    = nan(length(paras), 1);

for nGroup = 1:length(expression)
    
    indexExpression = strcmp(expression{nGroup}, {paras.expression});
    indexCaInd      = strcmp(CaIndicator{nGroup}, {paras.CaIndicator});
    
    group(indexExpression & indexCaInd)     = nGroup;
end


nTitles = {'\tau_{r}', '\tau_{d}', 'n', 'K', 'Fm'};

groupColor = [         0    0.4470    0.7410
    0.9290    0.6940    0.1250
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

figure;
[~, ax, ~] = gplotmatrix(paraMat, [], group, groupColor, 'oooo', [], 'off', [], nTitles, nTitles);

setPrint(8*2, 8*2, [PlotDir 'ModelCellFits/parasSpace'])