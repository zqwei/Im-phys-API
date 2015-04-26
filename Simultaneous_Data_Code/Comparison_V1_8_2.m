%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_8

addpath('../Func');
setDir;
load ([TempDatDir 'DataList.mat']);
addNoise         = [1 0 0 0 0 0];

if ~exist([PlotDir '/Collected_Units_Decodability_Epoch'],'dir')
    mkdir([PlotDir '/Collected_Units_Decodability_Epoch'])
end
numTrials           = 400;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = true(numTrainingTrials,1);
testTargets         = true(numTestTrials,1);
totTargets          = [testTargets; trainingTargets];

numFold               = 30;

for nData             = fileToAnalysis
    load([TempDatDir DataSetList(nData).name '.mat'])   
    figure;
    timePoints   = timePointTrialPeriod(DataSetList(nData).params.polein, DataSetList(nData).params.poleout, DataSetList(nData).params.timeSeries);
    numPeriods   = length(timePoints) - 1;
    decodability = zeros(numFold, numPeriods, size(nDataSet(1).unit_yes_trial,2));
    for nFold       = 1:numFold
        numUnits     = length(nDataSet);
        nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, numRandPickUnits)), totTargets, numTestTrials);
        % nSessionData : Ntrial x Nneuron x Nt
        nSessionData = permute(nSessionData,[1 3 2]);
        % nSessionData : Ntrial x Nt x Nneuron
        % numTrials    = size(nSessionData ,1);
        EpochIndex   = epochIndex(DataSetList(nData).params);
        EpochIndex   = EpochIndex(:,ones(1,numTrials))';
        % EpochIndex : Ntrial x Nt
        % numTestTrials = round(numTrials*0.2);        
        decodability(nFold,:,:)    = decodabilitySliceDataTaperAccumulated(nSessionData, EpochIndex, 0, addNoise(nData), numTestTrials, numPeriods);
    end
    hold on
    area(DataSetList(nData).params.timeSeries,squeeze(mean(decodability,1))','Edgecolor','none');
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0 1])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[0.5], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    ylabel('Decodability')
    xlabel('Time (s)')
    setPrint(4, 3, [PlotDir 'Collected_Units_Decodability_Epoch/Collected_Units_Decodability_EpochLDA1AllAccFixedNumUnitsYesTrial_' DataSetList(nData).name], 'pdf')
end


numTrials           = 400;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = false(numTrainingTrials,1);
testTargets         = false(numTestTrials,1);
totTargets          = [testTargets; trainingTargets];

numFold               = 30;

for nData             = fileToAnalysis
    load([TempDatDir DataSetList(nData).name '.mat'])   
    figure;
    timePoints   = timePointTrialPeriod(DataSetList(nData).params.polein, DataSetList(nData).params.poleout, DataSetList(nData).params.timeSeries);
    numPeriods   = length(timePoints) - 1;
    decodability = zeros(numFold, numPeriods, size(nDataSet(1).unit_yes_trial,2));
    for nFold       = 1:numFold
        numUnits     = length(nDataSet);
        nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, numRandPickUnits)), totTargets, numTestTrials);
        % nSessionData : Ntrial x Nneuron x Nt
        nSessionData = permute(nSessionData,[1 3 2]);
        % nSessionData : Ntrial x Nt x Nneuron
        % numTrials    = size(nSessionData ,1);
        EpochIndex   = epochIndex(DataSetList(nData).params);
        EpochIndex   = EpochIndex(:,ones(1,numTrials))';
        % EpochIndex : Ntrial x Nt
        % numTestTrials = round(numTrials*0.2);        
        decodability(nFold,:,:)    = decodabilitySliceDataTaperAccumulated(nSessionData, EpochIndex, 0, addNoise(nData), numTestTrials, numPeriods);
    end
    hold on
    area(DataSetList(nData).params.timeSeries,squeeze(mean(decodability,1))','Edgecolor','none');
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0 1])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[0.5], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    ylabel('Decodability')
    xlabel('Time (s)')
    setPrint(4, 3, [PlotDir 'Collected_Units_Decodability_Epoch/Collected_Units_Decodability_EpochLDA1AllAccFixedNumUnitsNoTrial_' DataSetList(nData).name], 'pdf')
end

close all
