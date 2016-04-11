%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decoding which trial period one is in, how fast can you tell the data
% that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using Shaul's layout
% https://www.evernote.com/l/AMwcR-KzK09KCrxHDHCMYjce1uXmXzHN8Co
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('../Func');
setDir;
load ([TempDatDir 'DataListEphys.mat']);
addNoise         = [1 1 1 1 1 1];

detectThres = 0.2;
cmap = [ 0    0.4470    0.7410
    0.6350    0.0780    0.1840
    0.4660    0.6740    0.1880
    0.8500    0.3250    0.0980
    0.4940    0.1840    0.5560
    0.3010    0.7450    0.9330];

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

nDatas          = 1:length(DataSetList);
meanDelays      = nan(length(nDatas), 3); % #nData, #Roc, #Epoch
semDelays       = nan(length(nDatas), 3);  
ROCThres        = 0.5;

for nData   = nDatas
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
    delayTimes = nan(numFold, 3);
    for nEpoch = 1:3
        preEpoch = nEpoch;
        newEpoch = nEpoch + 1;
        decodabilityEpoch      = decodability(:, :, timePoints(nEpoch+1):timePoints(nEpoch+2));
        delayTimes (:, nEpoch) = sum(decodabilityEpoch(:, nEpoch, :)>detectThres, 3);
    end
    delayTimes = delayTimes/DataSetList(nData).params.frameRate; 
    semDelay   = std(delayTimes)/sqrt(numFold);
    meanDelay  = mean(delayTimes, 1);
    semDelays(nData, :)  = semDelay*1000;
    meanDelays(nData, :) = meanDelay*1000;
end

figure;
labelDelays         = {'Sample', 'Delay', 'Response'};
for nPlot           = 1:length(labelDelays)
    subplot(1, length(labelDelays), nPlot)
    errorbar(nDatas, meanDelays(:, nPlot), semDelays(:, nPlot), 'ok')
    xlim([0.5 6.5])
    xlabel('Ephy boxcar filter (ms)')
    ylabel('latency to epoch change (ms)')
    set(gca, 'XTick', nDatas, 'XTickLabel', {'no 67', '70', '100', '150', '200', '250'})
    xticklabel_rotate([],30)
    set(gca, 'TickDir', 'out')
    title(labelDelays{nPlot});
    box off
end

setPrint(8*3, 6, [PlotDir 'EphysBinningTest/CollectedUnitsDecodabilityEpoch'])

