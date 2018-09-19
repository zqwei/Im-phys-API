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
ActiveNeuronIndex  = DataSetList(1).ActiveNeuronIndex';
params             = DataSetList(1).params;
DataSetListName{1} = DataSetList(1).name;
load ([TempDatDir 'DataListS2CModel.mat']);
DataSetListName{2} = DataSetList(3).name;
DataSetListName{3} = DataSetList(4).name;
addNoise           = [0 0 0 0];
timePoints         = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
EpochIndex         = epochIndex(params);

detectThres = 0.3;
cmap = [ 0    0.4470    0.7410
    0.6350    0.0780    0.1840
    0.4660    0.6740    0.1880];

numRandPickUnits    = 100;
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

ROCValue        = 0.5:0.10:0.80;
meanDelays      = nan(length(DataSetListName), length(ROCValue), 3); % #nData, #Roc, #Epoch
semDelays       = nan(length(DataSetListName), length(ROCValue), 3);  

for nROCThres   = 1:length(ROCValue)
    ROCThres    = ROCValue(nROCThres);
    for nData   = 1:length(DataSetListName)
        load([TempDatDir DataSetListName{nData} '.mat'])   
        selectedNeuronalIndex = ActiveNeuronIndex';
        selectedNeuronalIndex = selectedHighROCneurons(nDataSet, params, ROCThres, selectedNeuronalIndex);
        nDataSet              = nDataSet(selectedNeuronalIndex);
        numPeriods            = length(timePoints) - 1;
        decodability          = zeros(numFold, numPeriods, size(nDataSet(1).unit_yes_trial,2));
        for nFold        = 1:numFold
            numUnits     = length(nDataSet);
            nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, numRandPickUnits)), totTargets, numTestTrials);
            nSessionData = permute(nSessionData,[1 3 2]);            
            EpochIndexs  = EpochIndex(:,ones(1,numTrials))';
            decodability(nFold,:,:)    = decodabilitySliceDataTaperAccumulated(nSessionData, EpochIndexs, 0, addNoise(nData), numTestTrials, numPeriods);
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
        delayTimes = delayTimes/params.frameRate; 
        % epochTime  = [DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0];    
        semDelay   = std(delayTimes)/sqrt(numFold);
        % meanDelay  = mean(delayTimes, 1) - epochTime;
        meanDelay  = mean(delayTimes, 1);
        semDelays(nData, nROCThres, :)  = semDelay*1000;
        meanDelays(nData, nROCThres, :) = meanDelay*1000;
    end
end

figure;
labelDelays         = {'Sample', 'Delay', 'Response'};
for nPlot           = 1:length(labelDelays)
    subplot(1, length(labelDelays), nPlot)
    hold on
    for nData       = 1:length(DataSetListName)
        shadedErrorBar(ROCValue, squeeze(meanDelays(nData, :, nPlot)), ...
            squeeze(semDelays(nData, :, nPlot)), {'-','color',cmap(nData, :), 'linewid', 1.0}, 0.5)
    end
    xlim([0.45, 0.85])
%     ylim([-60 400])
    xlabel('ROC thres.')
    ylabel('latency to epoch change (ms)')
    set(gca, 'TickDir', 'out')
    title(labelDelays{nPlot});
    box off
end

subplot(1,3,1)
ylim([100 500])
subplot(1,3,2)
ylim([0 400])
subplot(1,3,3)
ylim([50 250])
setPrint(8*3, 6, [PlotDir 'S2CModel/CollectedUnitsDecodabilityEpochROC'])

numRandPickUnitsSet     = [50 100 200 500];% [20 50 100 200 500];
nDatas                  = 1:length(DataSetList);
meanDelays              = nan(length(DataSetListName), length(numRandPickUnitsSet), 3); % #nData, #Roc, #Epoch
semDelays               = nan(length(DataSetListName), length(numRandPickUnitsSet), 3);  

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
    for nData           = 1:length(DataSetListName)
        load([TempDatDir DataSetListName{nData} '.mat'])   
        selectedNeuronalIndex = ActiveNeuronIndex';
        selectedNeuronalIndex = selectedHighROCneurons(nDataSet, params, ROCThres, selectedNeuronalIndex);
        nDataSet              = nDataSet(selectedNeuronalIndex);
        numPeriods            = length(timePoints) - 1;
        decodability          = zeros(numFold, numPeriods, size(nDataSet(1).unit_yes_trial,2));
        for nFold        = 1:numFold
            numUnits     = length(nDataSet);
            nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, numRandPickUnits)), totTargets, numTestTrials);
            nSessionData = permute(nSessionData,[1 3 2]);
            EpochIndexs   = EpochIndex(:,ones(1,numTrials))';
            decodability(nFold,:,:)    = decodabilitySliceDataTaperAccumulated(nSessionData, EpochIndexs, 0, addNoise(nData), numTestTrials, numPeriods);
        end
%         [~, maxIndex] = max(decodability, [], 2);
%         maxIndex = squeeze(maxIndex);
%         delayTimes = nan(numFold, 3);
%         for nEpoch = 1:3
%             delayTimes (:, nEpoch) = sum(maxIndex < nEpoch+1, 2);
%         end
%         delayTimes = DataSetList(nData).params.timeSeries(delayTimes+1);    
        delayTimes = nan(numFold, 3);
        for nEpoch = 1:3
            preEpoch = nEpoch;
            newEpoch = nEpoch + 1;
            decodabilityEpoch      = decodability(:, :, timePoints(nEpoch+1):timePoints(nEpoch+2));
            delayTimes (:, nEpoch) = sum(decodabilityEpoch(:, nEpoch, :)>detectThres, 3);
        end
        delayTimes = delayTimes/params.frameRate; 
        semDelay   = std(delayTimes)/sqrt(numFold);
        meanDelay  = mean(delayTimes, 1);
        semDelays(nData, nRandPickUnits, :)  = semDelay*1000;
        meanDelays(nData, nRandPickUnits, :) = meanDelay*1000;
    end
end

figure;
labelDelays         = {'Sample', 'Delay', 'Response'};
for nPlot           = 1:length(labelDelays)
    subplot(1, length(labelDelays), nPlot)
    hold on
    for nData       = [1 2 3]
        shadedErrorBar(numRandPickUnitsSet, squeeze(meanDelays(nData, :, nPlot)), ...
            squeeze(semDelays(nData, :, nPlot)), {'-','color',cmap(nData, :), 'linewid', 1.0}, 0.5)
    end
    xlim([0, 501])
%     ylim([-60 400])
    xlabel('Num. units')
    ylabel('Latency to epoch change (ms)')
    set(gca, 'TickDir', 'out')
    title(labelDelays{nPlot});
    box off
end
subplot(1,3,1)
ylim([100 600])
subplot(1,3,2)
ylim([0 400])
subplot(1,3,3)
ylim([50 300])

setPrint(8*3, 6, [PlotDir 'S2CModel/CollectedUnitsDecodabilityEpochNumUnits'])