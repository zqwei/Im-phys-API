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

load ([TempDatDir 'DataListShuffle.mat']);
addNoise            = [1 0 0 0];

cmap                = cbrewer('div', 'Spectral', 128, 'cubic');


if ~exist([PlotDir '/CollectedUnitsDecodability'],'dir')
    mkdir([PlotDir '/CollectedUnitsDecodability'])
end

stepSize            = 10;

for nData             = [1 3 4]
    if nData == 1
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList     = false(length(nDataSet),1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end

    oldDataSet               = nDataSet;
    maxRandPickUnits         = 20;
    decodabilityAll          = zeros(maxRandPickUnits, numFold, size(nDataSet(1).unit_yes_trial,2));
    figure;
    hold on
    for numRandPickUnits      = 1:maxRandPickUnits;
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex(~neuronRemoveList)';% & [DataSetList(nData).cellinfo(:).cellType] == 1;
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
        decodabilityAll(numRandPickUnits, :, :) = decodability;
    end
    save(['decodabilityAll_' DataSetList(nData).name], 'decodabilityAll')
end