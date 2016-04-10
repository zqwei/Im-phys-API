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

ROCValue        = 0.5:0.1:0.8;
nDatas          = [1 3 4];
meanDelays      = nan(length(nDatas), length(ROCValue), 3); % #nData, #Roc, #Epoch
semDelays       = nan(length(nDatas), length(ROCValue), 3);  

for nROCThres   = 1:length(ROCValue)
    ROCThres    = ROCValue(nROCThres);
    for mData   = 1:length(nDatas)
        nData   = nDatas(mData);
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
%         [~, maxIndex] = max(decodability, [], 2);
%         maxIndex = squeeze(maxIndex);
        delayTimes = nan(numFold, 3);
        for nEpoch = 1:3
            delayTimes (:, nEpoch) = sum(maxIndex < nEpoch+1, 2);
        end
        delayTimes = DataSetList(nData).params.timeSeries(delayTimes+1);    
        epochTime  = [DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0];    
        semDelay   = std(delayTimes)/sqrt(numFold);
        meanDelay  = mean(delayTimes, 1) - epochTime;
        semDelays(mData, nROCThres, :)  = semDelay*1000;
        meanDelays(mData, nROCThres, :) = meanDelay*1000;
    end
end

figure;
labelDelays         = {'Stim.', 'Delay', 'Resp.'};
for nPlot           = 1:length(labelDelays)
    subplot(1, length(labelDelays), nPlot)
    hold on
    for nData       = 1:size(meanDelays, 1)
        shadedErrorBar(ROCValue, squeeze(meanDelays(nData, :, nPlot)), ...
            squeeze(semDelays(nData, :, nPlot)), {'-','color',cmap(nData, :), 'linewid', 1.0}, 0.5)
    end
    xlim([0.45, 0.85])
    ylim([-60 400])
    xlabel('ROC thres.')
    ylabel('latency to epoch change (ms)')
    set(gca, 'TickDir', 'out')
    title(labelDelays{nPlot});
end

setPrint(8*3, 6, [PlotDir 'CollectedUnitsDecodabilityEpoch/CollectedUnitsDecodabilityEpochROC'])

numRandPickUnitsSet     = [20 50 100 200 500];
nDatas                  = [1 3 4];
meanDelays              = nan(length(nDatas), length(numRandPickUnitsSet), 3); % #nData, #Roc, #Epoch
semDelays               = nan(length(nDatas), length(numRandPickUnitsSet), 3);  

for nRandPickUnits      = 1:length(numRandPickUnitsSet)
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
    for mData           = 1:length(nDatas)
        nData           = nDatas(mData);
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
        delayTimes = DataSetList(nData).params.timeSeries(delayTimes+1);    
        epochTime  = [DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0];    
        semDelay   = std(delayTimes)/sqrt(numFold);
        meanDelay  = mean(delayTimes, 1) - epochTime;
        semDelays(mData, nRandPickUnits, :)  = semDelay*1000;
        meanDelays(mData, nRandPickUnits, :) = meanDelay*1000;
    end
end

figure;
labelDelays         = {'Stim.', 'Delay', 'Resp.'};
for nPlot           = 1:length(labelDelays)
    subplot(1, length(labelDelays), nPlot)
    hold on
    for nData       = 1:size(meanDelays, 1)
        shadedErrorBar(numRandPickUnitsSet, squeeze(meanDelays(nData, :, nPlot)), ...
            squeeze(semDelays(nData, :, nPlot)), {'-','color',cmap(nData, :), 'linewid', 1.0}, 0.5)
    end
    xlim([0, 501])
    ylim([-60 400])
    xlabel('Num. units')
    ylabel('Latency to epoch change (ms)')
    set(gca, 'TickDir', 'out')
    title(labelDelays{nPlot});
end

setPrint(8*3, 6, [PlotDir 'CollectedUnitsDecodabilityEpoch/CollectedUnitsDecodabilityEpochNumUnits'])


margNames = {'Spike', 'GP4.3', '6s-AAV'};
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
