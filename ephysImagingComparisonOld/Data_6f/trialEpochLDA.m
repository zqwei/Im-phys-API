%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decoding which trial period one is in, how fast can you tell the data
% that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
addNoise         = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

if ~exist([PlotDir '/CollectedUnitsDecodabilityEpoch'],'dir')
    mkdir([PlotDir '/CollectedUnitsDecodabilityEpoch'])
end

cmap = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.4660    0.6740    0.1880
    0.6350    0.0780    0.1840];

numRandPickUnits    = 50;
numTrials           = numRandPickUnits*5;
numTestTrials       = numRandPickUnits*2;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
trainingTargets     = trainingTargets(randperm(numTrainingTrials));
testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
testTargets         = testTargets(randperm(numTestTrials));
totTargets          = [testTargets; trainingTargets];
ROCThres            = 0.5;
numFold             = 10;

for nData           = [1 10]
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat']);
        neuronRemoveList = false(length(nDataSet), 1);
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat']);
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex(~neuronRemoveList)';
    end
    
    depth_list          = [nDataSet.depth_in_um]';    
    oldDataSet          = nDataSet;
    selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    selectedNeuronalIndex = selectedNeuronalIndex & depth_list < 471;
    nDataSet              = oldDataSet(selectedNeuronalIndex);

    
    figure;
    timePoints            = timePointTrialPeriod(DataSetList(nData).params.polein, DataSetList(nData).params.poleout, DataSetList(nData).params.timeSeries);
    numPeriods            = length(timePoints) - 1;
    decodability          = zeros(numFold, numPeriods, size(nDataSet(1).unit_yes_trial,2));
    for nFold        = 1:numFold
        numUnits     = length(nDataSet);
        nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, numRandPickUnits)), totTargets, numTestTrials);
        nSessionData = permute(nSessionData,[1 3 2]);
        EpochIndex   = epochIndex(DataSetList(nData).params);
        EpochIndex   = EpochIndex(:,ones(1,numTrials))';
        decodability(nFold,:,:)    = decodabilitySliceDataTaperAccumulated(nSessionData, EpochIndex, 0, addNoise(nData), numTestTrials, numPeriods);
    end
   
    subplot(10, 1, 1)
    [~, maxIndex] = max(mean(decodability, 1), [], 2);
    maxIndex = squeeze(maxIndex);
    imagesc(1, DataSetList(nData).params.timeSeries, maxIndex', [1 4]);
    axis off
    colormap(cmap(1:4, :));
    
    subplot(10, 1, 2:9)
    hold on
    area(DataSetList(nData).params.timeSeries,squeeze(mean(decodability,1))','Edgecolor','none');
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0 1])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[ 0.5 ], 'Color','k','Linestyle','--','linewid', 0.5) %#ok<NBRAK>
    box off;
    hold off;
    ylabel('Decodability')
    xlabel('Time (s)')
    colormap(cmap(1:4, :));
    
    setPrint(8, 6, [PlotDir 'CollectedUnitsDecodabilityEpoch/CollectedUnitsDecodabilityEpoch_0_5_' DataSetList(nData).name])
end


numRandPickUnits = 10;
numTrials           = 500;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
trainingTargets     = trainingTargets(randperm(numTrainingTrials));
testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
testTargets         = testTargets(randperm(numTestTrials));
totTargets          = [testTargets; trainingTargets];
ROCThres            = 0.7;
numFold             = 10;

for nData           = 10
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex(~neuronRemoveList)';
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex(~neuronRemoveList)';
    end
    params      = DataSetList(nData).params;
    numTime     = length(params.timeSeries);
    oldDataSet          = nDataSet;
    depth               = [nDataSet.depth_in_um];
    selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);

    decodabilityDepth   = nan(8, numTime);

    for nDepth          = 0:100:700
        depthIndex      = depth>nDepth & depth<(nDepth + 100);        
        if sum(depthIndex)>10
            decodability    = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));    
            for nFold       = 1:numFold
                nDataSet            = oldDataSet(depthIndex);
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
            if numFold>1
                decodabilityDepth(nDepth/100+1, :) = mean(decodability);
            else
                decodabilityDepth(nDepth/100+1, :) = decodability;
            end
        end
    end
    
    figure;
    hold on;
    imagesc(DataSetList(nData).params.timeSeries, 0:-100:-700, decodabilityDepth, [0.5 1])
    axis xy;
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    set(gca, 'yTick', -800:400:0)
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('# units');
    colormap(cmap)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'ConfoundingFactorPCA/CollectedUnitsLDADepth_' DataSetList(nData).name])
end