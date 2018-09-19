addpath('../Func');
setDir;

load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
params            = DataSetList(1).params;

dataSetNames   = {'Modeled_6s_AAV_nRSParaSet_', 'Modeled_GP43_nRSParaSet_', 'Modeled_GP517_nRSParaSet_'};


for nData         = 1:3
    analysisMat   = repmat(struct('nParaSet',1, 'pCa0', 1, ...
                                'pTaud', 1, 'sizeGroup', 1),1000, 1);
    for nParaSet  = 1:1000
        load([TempDatDir 'Results_' dataSetNames{nData} num2str(nParaSet, '%04d') '.mat'], 'sizeGroup', 'countMaxId', 'peakiness', 'PCAVar', 'decodability')
        analysisMat(nParaSet).sizeGroup    = sizeGroup;
        analysisMat(nParaSet).peakiness    = peakiness;
        analysisMat(nParaSet).PCAVar       = PCAVar;
        analysisMat(nParaSet).decodability = decodability;
    end
    
    save([TempDatDir 'ResultsCompiled_' dataSetNames{nData} '.mat'], 'analysisMat')
end