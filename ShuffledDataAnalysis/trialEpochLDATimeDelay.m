%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decoding which trial period one is in, how fast can you tell the data
% that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
addNoise         = [1 0 0 0];

if ~exist([PlotDir '/CollectedUnitsDecodabilityEpoch'],'dir')
    mkdir([PlotDir '/CollectedUnitsDecodabilityEpoch'])
end

cmap = [ 0    0.4470    0.7410
    0.6350    0.0780    0.1840
    0.4660    0.6740    0.1880];

numRandPickUnits    = 100;
numTestTrials       = 200;
numFold             = 30;

numTrials           = numRandPickUnits * 3 + numTestTrials;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
trainingTargets     = trainingTargets(randperm(numTrainingTrials));
testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
testTargets         = testTargets(randperm(numTestTrials));
totTargets          = [testTargets; trainingTargets];


figure;


for ROCThres        = 0.5:0.1:0.8
    
    
    subplot(2, 4, ROCThres*10 - 4)
    
    meanDelays          = [];
    semDelays           = [];
    labelDelays         = {'Stim.', 'Delay', 'Resp.'};

    for nData           = [1 3 4]
        load([TempDatDir DataSetList(nData).name '.mat'])   
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
        selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
        nDataSet              = nDataSet(selectedNeuronalIndex);
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

        [~, maxIndex] = max(decodability, [], 2);
        maxIndex = squeeze(maxIndex);

        delayTimes = nan(numFold, 3);

        for nEpoch = 1:3
            delayTimes (:, nEpoch) = sum(maxIndex < nEpoch+1, 2);
        end

        hold on;
        delayTimes = DataSetList(nData).params.timeSeries(delayTimes+1);    
        epochTime  = [DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0];    
        semDelay   = std(delayTimes)/sqrt(numFold);
        meanDelay  = mean(delayTimes, 1) - epochTime;
        semDelays  = [semDelays, semDelay'*1000]; %#ok<AGROW>
        meanDelays = [meanDelays, meanDelay'*1000]; %#ok<AGROW>
    end

    barerror((1:3)', meanDelays, semDelays, 1, cmap(1:3, :), ['k', 'k', 'k'], labelDelays);
    ylabel('Delay to epoch (ms)')
    ylim([-60 400])
    xlim([0.5 3.5])
    set(gca, 'TickDir', 'out')
    title(['ROC =  ' num2str(ROCThres)]);
end

numRandPickUnitsSet     = [50 100 200 500];

for nRandPickUnits      = 1:4
    
    
    ROCThres            = 0.5;
    numRandPickUnits    = numRandPickUnitsSet(nRandPickUnits);
    numTestTrials       = 200;
    numFold             = 30;

    numTrials           = numRandPickUnits * 3 + numTestTrials;
    numTrainingTrials   = numTrials - numTestTrials;
    trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
    trainingTargets     = trainingTargets(randperm(numTrainingTrials));
    testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
    testTargets         = testTargets(randperm(numTestTrials));
    totTargets          = [testTargets; trainingTargets];
    
    
    subplot(2, 4, 4+nRandPickUnits)
    
    meanDelays          = [];
    semDelays           = [];
    labelDelays         = {'Stim.', 'Delay', 'Resp.'};

    for nData           = [1 3 4]
        load([TempDatDir DataSetList(nData).name '.mat'])   
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
        selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
        nDataSet              = nDataSet(selectedNeuronalIndex);
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

        [~, maxIndex] = max(decodability, [], 2);
        maxIndex = squeeze(maxIndex);

        delayTimes = nan(numFold, 3);

        for nEpoch = 1:3
            delayTimes (:, nEpoch) = sum(maxIndex < nEpoch+1, 2);
        end

        hold on;
        delayTimes = DataSetList(nData).params.timeSeries(delayTimes+1);    
        epochTime  = [DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0];    
        semDelay   = std(delayTimes)/sqrt(numFold);
        meanDelay  = mean(delayTimes, 1) - epochTime;
        semDelays  = [semDelays, semDelay'*1000]; %#ok<AGROW>
        meanDelays = [meanDelays, meanDelay'*1000]; %#ok<AGROW>
    end

    barerror((1:3)', meanDelays, semDelays, 1, cmap(1:3, :), ['k', 'k', 'k'], labelDelays);
    ylabel('Delay to epoch (ms)')
    ylim([-60 400])
    xlim([0.5 3.5])
    set(gca, 'TickDir', 'out')
    title(['# units =  ' num2str(numRandPickUnits)]);
end

setPrint(8*4, 6*2, [PlotDir 'CollectedUnitsDecodabilityEpoch/CollectedUnitsDecodabilityEpoch'])

margNames = {'Spike', 'Transgenic', 'AAV'};
figure;
hold on
for nColor = 1:length(margNames)
    plot(0, nColor, 's', 'color', cmap(nColor,:), 'MarkerFaceColor',cmap(nColor,:),'MarkerSize', 8)
    text(1, nColor, margNames{nColor})
end
xlim([0 10])
hold off
axis off
setPrint(3, 3, [PlotDir 'CollectedUnitsDecodabilityEpoch/CollectedUnitsDecodabilityEpochLabel'])

close all
