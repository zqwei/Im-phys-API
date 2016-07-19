%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Different ROC
% Same number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
ActiveNeuronIndex  = DataSetList(1).ActiveNeuronIndex';
params             = DataSetList(1).params;
DataSetListName{1} = DataSetList(1).name;
load ([TempDatDir 'DataListS2CModel.mat']);
DataSetListName{2} = DataSetList(3).name;
load ([TempDatDir 'DataListS2CGP43Model.mat']);
DataSetListName{3} = DataSetList(2).name;

numFold             = 30;
addNoise         = [1 0 0 0];

cmap = [         0    0.4470    0.7410
    0.4940    0.1840    0.5560
    0.9290    0.6940    0.1250];

numRandPickUnits    = 100;
numTrials           = numRandPickUnits*5;
numTestTrials       = numRandPickUnits*2;
numTrainingTrials   = numTrials - numTestTrials;
ROCThres            = 0.5;

figure;

for nData             = 1:3
    load([TempDatDir DataSetListName{nData} '.mat'])
    selectedNeuronalIndex = ActiveNeuronIndex';
    oldDataSet          = nDataSet;
    hold on
    selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, params, ROCThres, selectedNeuronalIndex);
    nDataSet              = oldDataSet(selectedNeuronalIndex);
    decodability          = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));        
    for nFold    = 1:numFold
        trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
        trainingTargets     = trainingTargets(randperm(numTrainingTrials));
        testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
        testTargets         = testTargets(randperm(numTestTrials));
        totTargets          = [testTargets; trainingTargets];

        trainingDecisions   = trainingTargets(randperm(numTrainingTrials));
        testDecisions       = testTargets(randperm(numTestTrials));
        totDecisions        = [testDecisions; trainingDecisions];

        randPickUnits       = randperm(length(nDataSet));
        randPickUnits       = randPickUnits(1:numRandPickUnits);

        nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
        decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
    end
    
    meandecodability = mean(decodability,1);
%     meanValue        = mean(meandecodability(1:8));
%     meandecodability = (meandecodability - meanValue)/(1-meanValue)*0.5+0.5;
    
    shadedErrorBar(params.timeSeries, meandecodability,...
        std(decodability, 1)/sqrt(numFold),...
        {'-', 'linewid', 1.0, 'color', cmap(nData,:)}, 0.5);  
    xlim([min(params.timeSeries) max(params.timeSeries)]);
    ylim([0.5 1])
    gridxy ([params.polein, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    set(gca, 'TickDir', 'out')
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('Decodability');
end

setPrint(8, 6, [PlotDir 'S2CGP43Model/CollectedUnitsDecodabilityROC_Summary'])
margNames = {'Spike', 'S2C 6s-AAV', 'S2C GP4.3'};


close all;