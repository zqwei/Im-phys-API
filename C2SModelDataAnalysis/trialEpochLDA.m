%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decoding which trial period one is in, how fast can you tell the data
% that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListC2SModel.mat']);
addNoise         = [1 0 0 0];

if ~exist([PlotDir 'C2SModel'],'dir')
    mkdir([PlotDir 'C2SModel'])
end

cmap = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.4660    0.6740    0.1880
    0.6350    0.0780    0.1840];

numRandPickUnits    = 100;
numTrials           = 500;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
trainingTargets     = trainingTargets(randperm(numTrainingTrials));
testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
testTargets         = testTargets(randperm(numTestTrials));
totTargets          = [testTargets; trainingTargets];
ROCThres            = 0.5;
numFold             = 30;

for nData           = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])   
    selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = nDataSet(selectedNeuronalIndex);
    figure;
    timePoints            = timePointTrialPeriod(DataSetList(nData).params.polein, DataSetList(nData).params.poleout, DataSetList(nData).params.timeSeries);
    numPeriods            = length(timePoints) - 1;
    decodability          = zeros(numFold, numPeriods, size(nDataSet(1).unit_yes_trial,2));
    for nFold        = 1:numFold
        numUnits     = length(nDataSet);
        nSessionData = shuffleSessionDataVec(nDataSet(randperm(numUnits, numRandPickUnits)), totTargets, numTestTrials);
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
    
    setPrint(8, 6, [PlotDir 'C2SModel/CollectedUnitsDecodabilityEpoch_0_5_' DataSetList(nData).name])
end


ROCThres            = 0.7;
numFold             = 30;

for nData           = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])   
    selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = nDataSet(selectedNeuronalIndex);
    figure;
    timePoints            = timePointTrialPeriod(DataSetList(nData).params.polein, DataSetList(nData).params.poleout, DataSetList(nData).params.timeSeries);
    numPeriods            = length(timePoints) - 1;
    decodability          = zeros(numFold, numPeriods, size(nDataSet(1).unit_yes_trial,2));
    for nFold        = 1:numFold
        numUnits     = length(nDataSet);
        nSessionData = shuffleSessionDataVec(nDataSet(randperm(numUnits, numRandPickUnits)), totTargets, numTestTrials);
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
    
    setPrint(8, 6, [PlotDir 'C2SModel/CollectedUnitsDecodabilityEpoch_0_7_' DataSetList(nData).name])
end

close all
