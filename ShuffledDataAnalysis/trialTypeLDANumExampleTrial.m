%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Same ROC
% Different number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('../Func');
setDir;

numFold             = 100;
numTrials           = 500;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
ROCThres            = 0.5;

load ([TempDatDir 'DataListShuffle.mat']);
addNoise            = [1 0 0 0];

cmap                = cbrewer('div', 'Spectral', 128, 'cubic');


if ~exist([PlotDir '/CollectedUnitsDecodability'],'dir')
    mkdir([PlotDir '/CollectedUnitsDecodability'])
end

stepSize            = [50];%[5 10 15 50 100 150];

for nData             = 4%[1 3 4]
    if nData == 1
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList     = false(length(nDataSet),1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end
    oldDataSet               = nDataSet;
    figure;
    hold on
    for numRandPickUnits      = 1:length(stepSize);        
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
            randPickUnits       = randPickUnits(1:stepSize(numRandPickUnits));

            nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
            decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
        end
        figure;
        hold on
        plot(DataSetList(nData).params.timeSeries, decodability(1:10,:), '-', 'linewid', 0.5, 'color', [0.5 0.5 0.5]);
        plot(DataSetList(nData).params.timeSeries, mean(decodability), '-k', 'linewid', 2.0);
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([0.5 1])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        xlabel('Time (s)');
        ylabel('# units');
        set(gca, 'TickDir', 'out')
        setPrint(8, 6, [PlotDir 'CollectedUnitsDecodability/CollectedUnitsDecodabilityFixedROCThres_0_5_Example' DataSetList(nData).name '_numNeuron_' num2str(stepSize(numRandPickUnits))])
    end    
end

