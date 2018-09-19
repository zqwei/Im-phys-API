%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Different ROC
% Same number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;

numFold             = 30;
load ([TempDatDir 'DataListShuffle.mat']);
addNoise            = 1;

cmap = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

numRandPickUnits    = 50;
numTrials           = numRandPickUnits*5;
numTestTrials       = numRandPickUnits*2;
numTrainingTrials   = numTrials - numTestTrials;
ROCThres            = 0.5;

nData               = 1;
load([TempDatDir DataSetList(nData).name '.mat'])
selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';

oldDataSet          = nDataSet;

figure
hold on
for frThres = [0 1 4 10] % spike count in this case
    load(['validMat_' num2str(frThres, '%02d')], 'validMat')
    decodability    = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));    
    for nFold       = 1:numFold
        nDataSet              = oldDataSet(validMat(:, nFold));
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
        decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise, trainingTargets, testTargets);
    end
    plot(DataSetList(nData).params.timeSeries, mean(decodability,1), 'linewid', 1, 'color', cmap(mod(frThres+1, 7), :));
    plot(DataSetList(nData).params.timeSeries, mean(decodability,1)-std(decodability, 1)/sqrt(numFold), 'linewid', 0.5, 'color', cmap(mod(frThres+1, 7), :));
    plot(DataSetList(nData).params.timeSeries, mean(decodability,1)+std(decodability, 1)/sqrt(numFold), 'linewid', 0.5, 'color', cmap(mod(frThres+1, 7), :));
end

load('caData.mat')
plot(DataSetList(nData).params.timeSeries, mean(caDecodability,1), 'linewid', 1, 'color', 'k');
plot(DataSetList(nData).params.timeSeries, mean(caDecodability,1)-std(caDecodability, 1)/sqrt(numFold), 'linewid', 0.5, 'color', 'k');
plot(DataSetList(nData).params.timeSeries, mean(caDecodability,1)+std(caDecodability, 1)/sqrt(numFold), 'linewid', 0.5, 'color', 'k');

xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
ylim([0.5 1])
gridxy ([DataSetList(nData).params.polein, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
set(gca, 'TickDir', 'out')
box off;
hold off;
xlabel('Time (s)');
ylabel('Decodability');
setPrint(8, 6, 'Decodability', 'pdf')

