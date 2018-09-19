%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Different ROC
% Same number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;

load ([TempDatDir 'DataListShuffle.mat']);
addNoise         = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

if ~exist([PlotDir '/CollectedUnitsDecodability'],'dir')
    mkdir([PlotDir '/CollectedUnitsDecodability'])
end

cmap                = cbrewer('div', 'Spectral', 128, 'cubic');

numFold             = 10;
numTrials           = 1000;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
ROCThres            = 0.70;
stepSize            = 10;

for nData             = [1 10]    %for 1 using maxRandPickUnits = 10, for 10 using maxRandPickUnits = 20
    if nData == 1
        load([TempDatDir DataSetList(nData).name '.mat'])
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex(~neuronRemoveList)';
    end
    oldDataSet          = nDataSet;
    selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    addTrialNumber        = mean(selectedNeuronalIndex)*4;
    selectedNeuronalIndex = double(selectedNeuronalIndex);
    selectedNeuronalIndex(selectedNeuronalIndex==0) = rand(sum(selectedNeuronalIndex==0), 1);
    selectedNeuronalIndex = selectedNeuronalIndex > 1 - addTrialNumber;
    depth_list            = [nDataSet.depth_in_um]';
    selectedNeuronalIndex = selectedNeuronalIndex & depth_list < 471;    
    nDataSet              = oldDataSet(selectedNeuronalIndex);
    maxRandPickUnits      = 20;
    decodabilityAll       = zeros(maxRandPickUnits, size(nDataSet(1).unit_yes_trial,2));    
    
    
    figure;
    hold on
    for numRandPickUnits      = 1:maxRandPickUnits
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
            decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials), trainingTargets, testTargets);
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
    colormap(cmap)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'CollectedUnitsDecodability/CollectedUnitsDecodabilityFixedROCThres_0_7_Fast_' DataSetList(nData).name])
end