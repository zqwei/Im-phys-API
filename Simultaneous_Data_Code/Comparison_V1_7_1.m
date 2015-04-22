%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.  Saturation of decodeability with number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_7

addpath('../Func');
setDir;
controlUnits        = 35;
numTrials           = 400;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
trainingTargets     = trainingTargets(randperm(numTrainingTrials));
testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
testTargets         = testTargets(randperm(numTestTrials));
totTargets          = [testTargets; trainingTargets];
load ([TempDatDir 'DataList.mat']);
addNoise            = [1 0 0 0 0 0];

if ~exist([PlotDir '/Collected_Units_LDA_Coeffs__KickOut'],'dir')
    mkdir([PlotDir '/Collected_Units_LDA_Coeffs__KickOut'])
end

perKickOut = 0:0.05:0.90;
numFold    = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.3  Plot of Sparsness using Kicking-out neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nData             = [1 3 4 5 6]%1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])  
    numUnits     = length(nDataSet);
    nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, controlUnits)), totTargets, numTestTrials);
    numUnits     = length(nDataSet);
    decodability = decodabilityLDAKickOut(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets, perKickOut, numFold);    
    figure;
    hold on
    imagesc(DataSetList(nData).params.timeSeries, perKickOut*100, decodability.mean);
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0 90])
    caxis([0.5 1])
    colorbar
    axis xy;
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('% KO neurons');
    setPrint(8, 6, [PlotDir 'Collected_Units_LDA_Coeffs__KickOut/Collected_Units_LDA_Coeffs_FixedNumUnits_KickOut_' DataSetList(nData).name], 'pdf')
end
