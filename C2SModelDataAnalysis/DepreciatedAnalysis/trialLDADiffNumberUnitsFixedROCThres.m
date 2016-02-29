%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Same ROC
% Different number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('../Func');
setDir;

numFold             = 10;
numTrials           = 1000;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
ROCThres            = 0.5;

load ([TempDatDir 'DataListModeled.mat']);
addNoise            = [1 0 0 0];

if ~exist([PlotDir '/ModeledCollectedUnitsDecodability'],'dir')
    mkdir([PlotDir '/ModeledCollectedUnitsDecodability'])
end

stepSize            = 10;

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    oldDataSet               = nDataSet;
    maxRandPickUnits         = 10;
    decodabilityAll          = zeros(maxRandPickUnits, size(nDataSet(1).unit_yes_trial,2));
    figure;
    hold on
    for numRandPickUnits      = 1:maxRandPickUnits;
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
        selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
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
            randPickUnits       = randPickUnits(1:numRandPickUnits*stepSize);

            nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
            decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
        end
        decodabilityAll(numRandPickUnits, :) = mean(decodability,1);
    end
    imagesc(DataSetList(nData).params.timeSeries, (1:maxRandPickUnits)*stepSize, decodabilityAll,[0.5 1]);
    axis xy;
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([stepSize maxRandPickUnits*stepSize])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('# units');
    colorbar
    setPrint(8, 6, [PlotDir 'ModeledCollectedUnitsDecodability/CollectedUnitsDecodabilityDiffNumberUnitsFixedROCThres_' DataSetList(nData).name], 'pdf')
end

close all;