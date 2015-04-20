%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9. Fraction of variance captured by PCA, per session, by trial period.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_9

addpath('../Func');
setDir;
load ([TempDatDir 'DataList.mat']);
addNoise         = [1 1 0 0 0];
numTrials           = 3000;
numTestTrials       = 600;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = rand(numTrainingTrials, 1) > 0.5;
testTargets         = rand(numTestTrials, 1) > 0.5;
totTargets          = [testTargets; trainingTargets];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10.1 Fraction of variance captured by PCA, per session, by trial period.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    numSession       = length(nDataSet3D);
    m                = ceil(sqrt(numSession + 1));
    figure;
    for nPlot        = 1:numSession
        subplot(m, m, nPlot)
        hold on
        nSessionData = [permute(nDataSet3D(nPlot).unit_yes_trial,[2, 1, 3]); permute(nDataSet3D(nPlot).unit_no_trial,[2, 1, 3])];
        nSessionData = normalizationDim(nSessionData, 2);
        perEV        = explainedPCA(nSessionData+randn(size(nSessionData))*1e-3/sqrt(numTrials)*addNoise(nData), 'CONSTANT', DataSetList(nData).params);
        plot(DataSetList(nData).params.timeSeries, perEV, 'k', 'linewid',1);
        perEV        = explainedPCA(nSessionData+randn(size(nSessionData))*1e-3/sqrt(numTrials)*addNoise(nData), 'TIME_VARYING', DataSetList(nData).params);
        plot(DataSetList(nData).params.timeSeries, perEV, 'r', 'linewid',1);
        perEV        = explainedPCA(nSessionData+randn(size(nSessionData))*1e-3/sqrt(numTrials)*addNoise(nData), 'STAGE_VARYING', DataSetList(nData).params);
        plot(DataSetList(nData).params.timeSeries, perEV, 'b', 'linewid',1);
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([0 1])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
    end

    subplot(m, m, nPlot+1)
    hold on
    nSessionData = shuffleSessionData(nDataSet, totTargets);
    nSessionData = normalizationDim(nSessionData, 2);
    perEV        = explainedPCA(nSessionData+randn(size(nSessionData))*1e-3/sqrt(numTrials)*addNoise(nData), 'CONSTANT', DataSetList(nData).params);
    plot(DataSetList(nData).params.timeSeries, perEV, 'k', 'linewid',1);
    perEV        = explainedPCA(nSessionData+randn(size(nSessionData))*1e-3/sqrt(numTrials)*addNoise(nData), 'TIME_VARYING', DataSetList(nData).params);
    plot(DataSetList(nData).params.timeSeries, perEV, 'r', 'linewid',1);
    perEV        = explainedPCA(nSessionData+randn(size(nSessionData))*1e-3/sqrt(numTrials)*addNoise(nData), 'STAGE_VARYING', DataSetList(nData).params);
    plot(DataSetList(nData).params.timeSeries, perEV, 'b', 'linewid',1);
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0 1])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    setPrint(m*4, m*3, [PlotDir 'Single_Session_PCA_' DataSetList(nData).name], 'pdf')
end