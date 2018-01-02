addpath('../Func');
setDir;

load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
params            = DataSetList(1).params;

dataSetNames   = {'Modeled_6s_AAV_nParaSet_', 'Modeled_GP43_nParaSet_'};


for nData         = 1:2
    analysisMat   = repmat(struct('nParaSet',1, 'pCa0', 1, ...
                                'pTaud', 1, 'sizeGroup', 1),900, 1);
    for nParaSet  = 1:900
        load([TempDatDir 'Results_' dataSetNames{nData} num2str(nParaSet, '%03d') '.mat'], 'sizeGroup', 'countMaxId', 'peakiness', 'PCAVar', 'decodability')
        analysisMat(nParaSet).nParaSet     = nParaSet;
        analysisMat(nParaSet).pCa0         = 10 + (mod(nParaSet-1, 30) + 1) * 3;
        analysisMat(nParaSet).pTaud        = 10 + floor((nParaSet-1)/30) * 3;
        analysisMat(nParaSet).sizeGroup    = sizeGroup;
        analysisMat(nParaSet).peakiness    = peakiness;
        analysisMat(nParaSet).PCAVar       = PCAVar;
        analysisMat(nParaSet).decodability = decodability;
    end
    
    save([TempDatDir 'ResultsCompiled_' dataSetNames{nData} '.mat'], 'analysisMat')
end