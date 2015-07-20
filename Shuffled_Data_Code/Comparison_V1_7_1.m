%
% Test for spiking dataset.
%
%
%


addpath('../Func');
setDir;

numRandPickUnits    = 200;
numTrials           = 5000;
numTestTrials       = 2000;
numTrainingTrials   = numTrials - numTestTrials;
% trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
% trainingTargets     = trainingTargets(randperm(numTrainingTrials));
% testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
% testTargets         = testTargets(randperm(numTestTrials));
% totTargets          = [testTargets; trainingTargets];
totTargets          = [true(numTrials/2,1)*2; false(numTrials/2,1)];
% totTargets(totTargets==true) = 2;
% totTargets(totTargets==false) = 1;
% 
% totTargets
% totTargetsChar      = ones(numTrials, 1);
% totTargetsChar(totTargets) = 'y';
% totTargetsChar(totTargets) = 'n';
ROCThres            = 0.7;
numFold             = 10;

minNumTrialToAnalysis  = 30;
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'None';
minFiringRate          = 10;
spikeDataSet           = getSpikeData(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize); %#ok<NASGU>

selectedNeuronalIndex  = findHighFiringUnits(spikeDataSet, params, minFiringRate);
selectedNeuronalIndex = selectedHighROCneurons(spikeDataSet, params, ROCThres, selectedNeuronalIndex);
nDataSet              = spikeDataSet(selectedNeuronalIndex);

numUnits     = length(nDataSet);
figure;
hold on
% nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, numRandPickUnits)), totTargets, numTrials);
nSessionData = shuffleSessionData(nDataSet, totTargets, numTrials);
% nSessionData = nSessionData + randn(size(nSessionData))*1e-6/sqrt(numTrials)* addNoise(nData);
% nSessionData = normalizationDim(nSessionData, 2);
% coeffs       = coeffLDA([nSessionData; nSessionData], totTargets, totTargets);
coeffs       = coeffLDA(nSessionData, totTargets);
imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, coeffs'*coeffs);
xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
caxis([-1 1]);
axis xy;
gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
box off;
hold off;
xlabel('Time (s)')
ylabel('Time (s)')

