%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decoding which trial period one is in, how fast can you tell the data
% that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using Shaul's layout
% https://www.evernote.com/l/AMwcR-KzK09KCrxHDHCMYjce1uXmXzHN8Co
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
addNoise         = [1 0 0 0];

if ~exist([PlotDir '/CollectedUnitsDecodabilityEpoch'],'dir')
    mkdir([PlotDir '/CollectedUnitsDecodabilityEpoch'])
end

detectThres = 0.2;
cmap = [ 0    0.4470    0.7410
    0.6350    0.0780    0.1840
    0.4660    0.6740    0.1880];

numRandPickUnits    = 50;
numTestTrials       = 200;
numFold             = 30;
numTrials           = numRandPickUnits * 3 + numTestTrials;
numTrials           = max(numTrials, numTestTrials*5);
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
trainingTargets     = trainingTargets(randperm(numTrainingTrials));
testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
testTargets         = testTargets(randperm(numTestTrials));
totTargets          = [testTargets; trainingTargets];

ROCValue        = 0.5:0.05:0.65;
nDatas          = [1 3 4];
meanDelays      = nan(length(nDatas), length(ROCValue), 3); % #nData, #Roc, #Epoch
semDelays       = nan(length(nDatas), length(ROCValue), 3);  

for nROCThres   = 1:length(ROCValue)
    ROCThres    = ROCValue(nROCThres);
    for mData   = 1:length(nDatas)
        nData   = nDatas(mData);
        load([TempDatDir DataSetList(nData).name '.mat'])   
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
        selectedNeuronalIndex = selectedHighLocalROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
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
        % [~, maxIndex] = max(decodability, [], 2);
        % maxIndex = squeeze(maxIndex);
        delayTimes = nan(numFold, 3);
        for nEpoch = 1:3
            preEpoch = nEpoch;
            newEpoch = nEpoch + 1;
            decodabilityEpoch      = decodability(:, :, timePoints(nEpoch+1):timePoints(nEpoch+2));
            delayTimes (:, nEpoch) = sum(decodabilityEpoch(:, nEpoch, :)>detectThres, 3);
            % delayTimes (:, nEpoch) = sum(maxIndex < nEpoch+1, 2);
        end
        % delayTimes = DataSetList(nData).params.timeSeries(delayTimes+1); 
        delayTimes = delayTimes/DataSetList(nData).params.frameRate; 
        % epochTime  = [DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0];    
        semDelay   = std(delayTimes)/sqrt(numFold);
        % meanDelay  = mean(delayTimes, 1) - epochTime;
        meanDelay  = mean(delayTimes, 1);
        semDelays(mData, nROCThres, :)  = semDelay*1000;
        meanDelays(mData, nROCThres, :) = meanDelay*1000;
    end
end

figure;
labelDelays         = {'Sample', 'Delay', 'Response'};
for nPlot           = 1:length(labelDelays)
    subplot(1, length(labelDelays), nPlot)
    hold on
    for nData       = 1:size(meanDelays, 1)
        shadedErrorBar(ROCValue, squeeze(meanDelays(nData, :, nPlot)), ...
            squeeze(semDelays(nData, :, nPlot)), {'-','color',cmap(nData, :), 'linewid', 1.0}, 0.5)
    end
    xlim([0.45, 0.75])
%     ylim([-60 400])
    xlabel('ROC thres.')
    ylabel('latency to epoch change (ms)')
    set(gca, 'TickDir', 'out')
    title(labelDelays{nPlot});
    box off
end

setPrint(8*3, 6, [PlotDir 'CollectedUnitsDecodabilityEpoch/CollectedUnitsDecodabilityEpochLocalROC50CellSession'])


numRandPickUnits    = 20;
numTestTrials       = 200;
numFold             = 30;
numTrials           = numRandPickUnits * 3 + numTestTrials;
numTrials           = max(numTrials, numTestTrials*5);
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
trainingTargets     = trainingTargets(randperm(numTrainingTrials));
testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
testTargets         = testTargets(randperm(numTestTrials));
totTargets          = [testTargets; trainingTargets];

ROCValue        = 0.5:0.05:0.70;
nDatas          = [1 3 4];
meanDelays      = nan(length(nDatas), length(ROCValue), 3); % #nData, #Roc, #Epoch
semDelays       = nan(length(nDatas), length(ROCValue), 3);  

for nROCThres   = 1:length(ROCValue)
    ROCThres    = ROCValue(nROCThres);
    for mData   = 1:length(nDatas)
        nData   = nDatas(mData);
        load([TempDatDir DataSetList(nData).name '.mat'])   
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
        selectedNeuronalIndex = selectedHighLocalROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
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
        % [~, maxIndex] = max(decodability, [], 2);
        % maxIndex = squeeze(maxIndex);
        delayTimes = nan(numFold, 3);
        for nEpoch = 1:3
            preEpoch = nEpoch;
            newEpoch = nEpoch + 1;
            decodabilityEpoch      = decodability(:, :, timePoints(nEpoch+1):timePoints(nEpoch+2));
            delayTimes (:, nEpoch) = sum(decodabilityEpoch(:, nEpoch, :)>detectThres, 3);
            % delayTimes (:, nEpoch) = sum(maxIndex < nEpoch+1, 2);
        end
        % delayTimes = DataSetList(nData).params.timeSeries(delayTimes+1); 
        delayTimes = delayTimes/DataSetList(nData).params.frameRate; 
        % epochTime  = [DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0];    
        semDelay   = std(delayTimes)/sqrt(numFold);
        % meanDelay  = mean(delayTimes, 1) - epochTime;
        meanDelay  = mean(delayTimes, 1);
        semDelays(mData, nROCThres, :)  = semDelay*1000;
        meanDelays(mData, nROCThres, :) = meanDelay*1000;
    end
end

figure;
labelDelays         = {'Sample', 'Delay', 'Response'};
for nPlot           = 1:length(labelDelays)
    subplot(1, length(labelDelays), nPlot)
    hold on
    for nData       = 1:size(meanDelays, 1)
        shadedErrorBar(ROCValue, squeeze(meanDelays(nData, :, nPlot)), ...
            squeeze(semDelays(nData, :, nPlot)), {'-','color',cmap(nData, :), 'linewid', 1.0}, 0.5)
    end
    xlim([0.45, 0.75])
%     ylim([-60 400])
    xlabel('ROC thres.')
    ylabel('latency to epoch change (ms)')
    set(gca, 'TickDir', 'out')
    title(labelDelays{nPlot});
    box off
end

setPrint(8*3, 6, [PlotDir 'CollectedUnitsDecodabilityEpoch/CollectedUnitsDecodabilityEpochLocalROC20CellSession'])

