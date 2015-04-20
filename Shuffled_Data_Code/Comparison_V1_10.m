%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Fraction of variance captured by PCA, per session, by trial period.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_10

function Comparison_V1_10

    numTrials           = 3000;
    numTestTrials       = 600;
    numTrainingTrials   = numTrials - numTestTrials;
    trainingTargets     = rand(numTrainingTrials, 1) > 0.5;
    testTargets         = rand(numTestTrials, 1) > 0.5;
    totTargets          = [testTargets; trainingTargets];
    load ('TempDat/DataList.mat');
    addNoise            = [1 1 0 0];        
    Comparison_V1_10_1(DataSetList, numTestTrials, totTargets, addNoise, numTrials)
    % Comparison_V1_10_2(DataSetList, numTestTrials, totTargets, addNoise, numTrials)
    % Comparison_V1_10_3(DataSetList, numTestTrials, totTargets, addNoise, numTrials)
    % Comparison_V1_10_4(DataSetList, numTestTrials, totTargets, addNoise, numTrials)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10.1 Fraction of variance captured by PCA, per session, by trial period.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Comparison_V1_10_1(DataSetList, numTestTrials, totTargets, addNoise, numTrials)    
    for nData             = 1:length(DataSetList)
        load(['TempDat/' DataSetList(nData).name '.mat'])
        
        sessionIndex     = [nDataSet(:).sessionIndex];
        sessionVec       = unique(sessionIndex);
        numSession       = length(sessionVec);
        m                = ceil(sqrt(numSession + 1));
        figure;
        for nPlot        = 1:numSession
            if sum (sessionIndex == sessionVec(nPlot)) > 1
                subplot(m, m, nPlot)
                hold on
                nSessionData = shuffleSessionData(nDataSet(sessionIndex == sessionVec(nPlot)), totTargets);
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
        end
        
        subplot(m, m, nPlot+1)
        hold on
        nSessionData = shuffleSessionData(nDataSet, totTargets);
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
                
        setPrint(m*4, m*3, ['Plot/Single_Session_PCA_' DataSetList(nData).name], 'pdf')
        
    end
end