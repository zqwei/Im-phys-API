addpath('../Func');
setDir;

load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
params        = DataSetList(3).params;

dataSetNames  = {DataSetList(4).name, DataSetList(3).name};

% performance in raw imaging data
% col 1: monophasic
% col 2: multiphasic
% col 3: peakiness
% col 4: pca
% col 5: lda

performanceMat(1).mono  = 0.66;
performanceMat(1).multi = 0.05;
performanceMat(1).peak  = 0.49;
performanceMat(2).pca   = [0.280973899350333;0.0522282453019869;0.113237355356368];
performanceMat(2).ldaS  = [0.600526315789474;0.604736842105263;0.584736842105263;0.658947368421053;0.514736842105263;0.609473684210526;0.652105263157895;0.619473684210526;0.574210526315790;0.692631578947369;0.607894736842105;0.620000000000000;0.541578947368421;0.611578947368421;0.605789473684210;0.624210526315789;0.608947368421053;0.638421052631579;0.727368421052632;0.701052631578947;0.648421052631579;0.613157894736842;0.608947368421053;0.710000000000000;0.662105263157895;0.566842105263158;0.570000000000000;0.515789473684211;0.676842105263158;0.598421052631579];
performanceMat(1).ldaD  = [0.926818181818182;0.843181818181818;0.825000000000000;0.915909090909091;0.712727272727273;0.935454545454545;0.924545454545455;0.862727272727273;0.928636363636364;0.938181818181818;0.960000000000000;0.872272727272727;0.901818181818182;0.812272727272727;0.879090909090909;0.842272727272727;0.917727272727273;0.844090909090909;0.989090909090909;0.946363636363636;0.904090909090909;0.879090909090909;0.919545454545455;0.960000000000000;0.971818181818182;0.824090909090909;0.824545454545455;0.889545454545455;0.960909090909091;0.900454545454546];

performanceMat(2).mono  = 0.50;
performanceMat(2).multi = 0.04;
performanceMat(2).peak  = 0.38;
performanceMat(1).pca   = [0.146199627211641;0.0643145311901328;0.100584045801060];
performanceMat(1).ldaS  = [0.619473684210526;0.637894736842105;0.581578947368421;0.730526315789474;0.692631578947368;0.478421052631579;0.664736842105263;0.744736842105263;0.648421052631579;0.761578947368421;0.597368421052631;0.732105263157895;0.615789473684211;0.651052631578947;0.644210526315790;0.676842105263158;0.656842105263158;0.546842105263158;0.725263157894737;0.651052631578947;0.598421052631579;0.694210526315789;0.628421052631579;0.596315789473684;0.554736842105263;0.648947368421053;0.636315789473684;0.748947368421053;0.698947368421053;0.586842105263158];
performanceMat(2).ldaD  = [0.855909090909091;0.946818181818182;0.826363636363636;0.974090909090909;0.945909090909091;0.938181818181818;0.817272727272728;0.921818181818182;0.941363636363636;0.945909090909091;0.903636363636364;0.949090909090909;0.844090909090909;0.821363636363637;0.939090909090909;0.813181818181818;0.758636363636364;0.780909090909091;0.874545454545455;0.884090909090909;0.839090909090909;0.857272727272727;0.920000000000000;0.764090909090909;0.846818181818182;0.669090909090909;0.691363636363636;0.921363636363636;0.970909090909091;0.882272727272727];

refEphys.pca            = [0.00180378089102542,0.224695451595608,0.00254116035619819];
refEphys.ldaS           = [0.696315789473684;0.621578947368421;0.678947368421053;0.648421052631579;0.717368421052632;0.663157894736842;0.599473684210527;0.677894736842106;0.687894736842105;0.660000000000000;0.645789473684211;0.693157894736842;0.659473684210526;0.612631578947368;0.656315789473684;0.727368421052632;0.674736842105263;0.643684210526316;0.658947368421053;0.663157894736842;0.627368421052632;0.693157894736842;0.626842105263158;0.665789473684211;0.724736842105263;0.600526315789474;0.627368421052632;0.691052631578947;0.626842105263158;0.722631578947369];
refEphys.ldaD           = [0.862272727272728;0.785000000000000;0.856363636363637;0.782272727272727;0.845454545454546;0.852272727272727;0.758181818181818;0.851818181818182;0.855909090909091;0.829090909090909;0.804545454545455;0.780909090909091;0.799545454545455;0.845454545454546;0.812727272727273;0.831818181818182;0.862272727272727;0.817272727272727;0.838636363636364;0.774090909090909;0.802727272727273;0.864545454545455;0.781363636363636;0.813636363636364;0.853181818181818;0.659545454545455;0.895454545454545;0.915000000000000;0.834090909090909;0.864545454545455];

numComps = 3;


for nData = 1:2
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    
    % selectivity
    Z         = reshape([analysisMat.sizeGroup], 3, 1000);
    numNeuron = sum(Z(:,1));
    disp(mean(Z')/numNeuron)
    disp(std(Z')/numNeuron)
    figure;
    subplot(1, 6, 1)
    hold on
    hist(Z(2,:)/numNeuron, 100);
    gridxy([performanceMat(nData).mono], [],'Color','r','Linestyle','--');
    gridxy([0.58], [],'Color','k','Linestyle','--');
    hold off
    box off
    xlim([0 0.7])
    xlabel('Fraction monophasic neuron')
    ylabel('# Units')
    set(gca, 'TickDir', 'out')

    subplot(1, 6, 2)
    hold on
    hist(Z(3,:)/numNeuron, 100);
    gridxy([performanceMat(nData).multi], [],'Color','r','Linestyle','--');
    gridxy([0.31], [],'Color','k','Linestyle','--');
    hold off
    box off
    xlim([0 0.35])
    xlabel('Fraction multiphasic neuron')
    ylabel('# Units')
    set(gca, 'TickDir', 'out')
    
    % peakiness
    Z     = [analysisMat.peakiness];
    disp(mean(Z'))
    disp(std(Z'))
    subplot(1, 6, 3)
    hist(Z, 100);
    gridxy([performanceMat(nData).peak], [],'Color','r','Linestyle','--');
    gridxy([1.27], [],'Color','k','Linestyle','--');
    hold off
    box off
    xlim([0 1.3])
    xlabel('Peakiness')
    ylabel('# Units')
    set(gca, 'TickDir', 'out')
    
    % pca
    pcaFracTrial  = nan(1000, 1);
    for   nParams = 1:1000
        PCAVar    = analysisMat(nParams).PCAVar;
        pcaFracTrial(nParams) = 1 - PCAVar(1, 2)/sum(PCAVar(1, :));
    end
    refPCAImage = 1 - performanceMat(nData).pca(2)/sum(performanceMat(nData).pca);
    refPCAEphys = 1 - refEphys.pca(2)/sum(refEphys.pca);    
    subplot(1, 6, 4)
    disp(mean(pcaFracTrial))
    disp(std(pcaFracTrial))
    hist(pcaFracTrial, 100);
    gridxy(refPCAImage, [],'Color','r','Linestyle','--');
    gridxy(refPCAEphys, [],'Color','k','Linestyle','--');
    hold off
    box off
    xlim([0 1])
    xlabel('Fraction of Non-time component in 1st PC')
    ylabel('# Units')
    set(gca, 'TickDir', 'out')
    
    
    % LDA Sample
    ldaDecodes    = nan(1000, 1);
    for   nParams = 1:1000
        ldaDecodes(nParams) = mean(mean(analysisMat(nParams).decodability(:, 8:26)));
    end
    subplot(1, 6, 5)
    disp(mean(ldaDecodes))
    disp(std(ldaDecodes))
    hist(ldaDecodes, 100);
    gridxy(mean(performanceMat(nData).ldaS), [],'Color','r','Linestyle','--');
    gridxy(mean(refEphys.ldaS), [],'Color','k','Linestyle','--');
    hold off
    box off
    xlim([0.5 1])
    xlabel('Sample decodability')
    ylabel('# Units')
    set(gca, 'TickDir', 'out')
    
    % LDA Delay
    ldaDecodes    = nan(1000, 1);
    for   nParams = 1:1000
        ldaDecodes(nParams) = mean(mean(analysisMat(nParams).decodability(:, 26:47)));
    end
    subplot(1, 6, 6)
    disp(mean(ldaDecodes))
    disp(std(ldaDecodes))
    hist(ldaDecodes, 100);
    gridxy(mean(performanceMat(nData).ldaD), [],'Color','r','Linestyle','--');
    gridxy(mean(refEphys.ldaD), [],'Color','k','Linestyle','--');
    hold off
    box off
    xlim([0.5 1])
    xlabel('Sample decodability')
    ylabel('# Units')
    set(gca, 'TickDir', 'out')
    
    
end