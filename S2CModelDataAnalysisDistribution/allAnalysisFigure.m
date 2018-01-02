addpath('../Func');
setDir;

load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
params            = DataSetList(1).params;

dataSetNames   = {'Modeled_6s_AAV_nParaSet_', 'Modeled_GP43_nParaSet_'};


for nData         = 1:2
    load([TempDatDir 'ResultsCompiled_' dataSetNames{nData} '.mat'], 'analysisMat')
    
    X             = reshape([analysisMat.pCa0], 30, 30);
    Y             = reshape([analysisMat.pTaud], 30, 30);
    Z             = reshape([analysisMat.sizeGroup], 3, 900);
    figure;
    subplot(1, 2, 1)
    pcolor(X, Y, reshape(Z(2,:), 30, 30)/720);
    caxis([0 1])
    subplot(1, 2, 2)
    pcolor(X, Y, reshape(Z(3,:), 30, 30)/720);
    caxis([0 1])
end