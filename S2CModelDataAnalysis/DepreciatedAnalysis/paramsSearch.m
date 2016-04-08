% 
% Comparison based on single unit acitivity
% Generating Ca++ imaging data from ephys data using Tsai-Wen's model
% 
% -------------------------------------------------------------------------
% version 1.0
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;

minNumTrialToAnalysis  = 20;
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
minFiringRate          = 5;
nDataSet               = getSpikeDataWithEphysTime(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);
nData_ActiveNeuronIndex = findHighFiringUnits(nDataSet, params, minFiringRate);

% Modeled GCaMP6s
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
params.K               = random(truncatedNormal, length(nDataSet), 1) *  5.8337 + 13.9248;
params.n               = random(truncatedNormal, length(nDataSet), 1) *  0.4106 +  1.7531;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  0.0344 +  0.0728;
params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  0.5055 +  1.4551;
params.Fm              = params.Fm;
params.tau_r           = params.tau_r;
params.tau_d           = params.tau_d*2;
params.intNoise        = 1.5;
params.extNoise        = 1.5;
nDataSet               = getFakeCaImagingData(nDataSet, params);

nData_name             = 'Modeled_Ca_01';
nData_params  = params; 
save([TempDatDir nData_name '.mat'], 'nDataSet', 'nData_params'); 
load([TempDatDir nData_name '.mat'], 'nDataSet', 'nData_params'); 

%%%%rawActivity
ylabels                 = 'dF/F';
yAxes_set               = [-0.5 2.0];
lowFiringThres          = 0.3;

plotMeanActivityImagescWithSortWithCellinfo(nDataSet, nData_params, [], [], ylabels, lowFiringThres, yAxes_set); 
setPrint(6*4, 3*3, [PlotDir 'ModeledSingleUnitsImagescWithSort/SingleUnitsImagescWithSort_' nData_name], 'tif')

%%%%rocEpochDist
plotROCPopAccLines(nDataSet, nData_params); 
setPrint(8, 6, [PlotDir 'ModeledSingleUnitsROC/SingleUnitsROC_' nData_name], 'pdf')

%%%%switchSelectivityDist
logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, nData_params);
save([TempDatDir 'LogPValueTscore_' nData_name '.mat'], 'logPValueEpoch')
load([TempDatDir 'LogPValueTscore_' nData_name '.mat'], 'logPValue', 'logPValueEpoch')
unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);
sizeGroup = histcounts(unitGroup, 0:3);
figure;
groupNames      = {'Non-selective', 'Homogenous', 'Heterogenous'};
pie(sizeGroup, groupNames)
title('Distribution of cell types')
setPrint(6, 4.5, [PlotDir 'ModeledSingleUnitsTscore/SingleUnitsTscore_' nData_name], 'pdf')

%%%%FPVT
% Gaussian filter for spiking data
sigma                         = 0.1 / nData_params.binsize; % 200 ms
filterLength                  = 10;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);
numUnits                  = length(nDataSet);
pValue                    = applyFuncToCompareTrialType(nDataSet, @pValueTTest2);
meanDiffValue             = applyFuncToCompareTrialType(nDataSet, @meanDiff);
logPValue                 = -log(pValue);
% yes  -- blue trial
% no   -- red trial
zScores                   = -sign(meanDiffValue).*logPValue;     
actMat                    = logPValue;
numT                      = size(zScores,2);
bumpActThres              = 3; % > bumpActThres considering as a bump % 3 = -log(0.05)
bumpMat                   = actMat > bumpActThres;
bumpSize                  = ones(size(actMat,1),1);
bumpStartPoint            = nan(size(actMat,1),1);
bumpSign                  = ones(size(actMat,1),1);

for nUnit = 1: size(actMat,1)
    diffActMat                = diff([bumpMat(nUnit,:),0]);
    beginPoint                = find(diffActMat==1);
    endPoint                  = find(diffActMat==-1);
    if length(endPoint)>length(beginPoint); endPoint = endPoint(2:end); end
    [bumplength, bumpIndex]   = max(endPoint - beginPoint);
    if isempty(bumplength); bumplength = 0; end
    bumpSize(nUnit)           = bumplength;
    if bumplength==0
        bumpStartPoint(nUnit) = numT;
        bumpSign(nUnit)       = -1;
    else 
        bumpStartPoint(nUnit) = beginPoint(bumpIndex);
        bumpSign(nUnit)       = actMat(nUnit, bumpStartPoint(nUnit))/zScores(nUnit, bumpStartPoint(nUnit));
    end
end
[~, similaritySort]           = sortrows([bumpStartPoint, bumpSize, bumpSign], [-3 -1 -2]);    

h = figure;
hold on
imagesc(nData_params.timeSeries, 1:numUnits, zScores(similaritySort,:));
gridxy ([nData_params.polein, nData_params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
xlim([nData_params.timeSeries(1) nData_params.timeSeries(end)])
ylim([1 numUnits])
colormap(h, french(128,2))
%     colorbar;
caxis([-5 5])
axis xy
xlabel('Time (s)')
ylabel('Neuron Index')
box off;
setPrint(8, 6, [PlotDir 'ModeledSingleUnitsFPVT/SingleUnitsZScore_' nData_name], 'pdf')


%%%%meanDPCA
combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numComps             = 15;
time               = nData_params.timeSeries;
timeEvents         = [nData_params.polein, nData_params.poleout, 0];
firingRates        = generateDPCAData(nDataSet, numTrials);
firingRatesAverage = nanmean(firingRates, ndims(firingRates));
trialNum           = ones(size(firingRatesAverage, 1), size(firingRatesAverage, 2))*numTrials;
optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, ...
            trialNum, ...
            'numComps', numComps, ....
            'combinedParams', combinedParams, ...
            'numRep', 10, ...  % increase this number to ~10 for better accuracy
            'filename', [TempDatDir 'optimalLambdas_' nData_name '.mat'],...
            'display','no');
[W,V,whichMarg] = dpca(firingRatesAverage, numComps, ...
            'combinedParams', combinedParams, ...
            'lambda', optimalLambda);
explVar = dpcaEV(firingRatesAverage, W, 'combinedParams', combinedParams);

save([TempDatDir 'dPCA_' nData_name '.mat'], 'optimalLambda', 'W', 'V', 'explVar')

load([TempDatDir 'dPCA_' nData_name '.mat'], 'explVar')

figure;
subplot(1, 3, 1)
plot(1:numComps, [explVar.PCAcomponentVar(1:numComps)', explVar.dPCAcomponentVar'], '-o');
legend({'PCA', 'dPCA'})
box off
xlim([0 numComps+1])
set(gca,'xTick', 0:5:numComps)
xlabel('Component index')
ylabel('% EV')

subplot(1, 3, 2)
bar(1:numComps, explVar.PCAmargVar(:, 1:numComps)', 'stack', 'edgecolor', 'none');
box off
legend(margNames)
xlim([0 numComps+1])
xlabel('Component index')
set(gca,'xTick', 0:5:numComps)
ylabel('% EV per PC')

subplot(1, 3, 3)
bar(1:numComps, explVar.dPCAmargVar', 'stack', 'edgecolor', 'none');
box off
legend(margNames)
xlim([0 numComps+1])
xlabel('Component index')
set(gca,'xTick', 0:5:numComps)
ylabel('% EV per dPC')

setPrint(8*3, 6, [PlotDir 'ModeledCollectedUnitsdPCA/CollectedUnitsdPCA_' nData_name], 'pdf')

figure
bar(1:numComps, explVar.PCAmargVar(:, 1:numComps)', 'stack', 'edgecolor', 'none');
box off
legend(margNames)
xlim([0 numComps+1])
xlabel('Component index')
set(gca,'xTick', 0:5:numComps)
ylabel('% EV per PC')
setPrint(8, 6, [PlotDir 'ModeledCollectedUnitsdPCA/CollectedUnitsPCA_' nData_name], 'pdf')

%%%%trialLDADiffNumberUnitsFixedROCThres
numFold             = 10;
numTrials           = 1000;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
ROCThres            = 0.5;
stepSize            = 10;
oldDataSet               = nDataSet;
maxRandPickUnits         = 10;
decodabilityAll          = zeros(maxRandPickUnits, size(nDataSet(1).unit_yes_trial,2));
figure;
hold on
for numRandPickUnits      = 1:maxRandPickUnits;
    selectedNeuronalIndex = nData_ActiveNeuronIndex';
    selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, nData_params, ROCThres, selectedNeuronalIndex);
    nDataSet              = oldDataSet(selectedNeuronalIndex);
    decodability          = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));

    for nFold    = 1:numFold
        trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
        trainingTargets     = trainingTargets(randperm(numTrainingTrials));
        testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
        testTargets         = testTargets(randperm(numTestTrials));
        totTargets          = [testTargets; trainingTargets];

        trainingDecisions   = trainingTargets(randperm(numTrainingTrials));
        testDecisions       = testTargets(randperm(numTestTrials));
        totDecisions        = [testDecisions; trainingDecisions];

        randPickUnits       = randperm(length(nDataSet));
        randPickUnits       = randPickUnits(1:numRandPickUnits*stepSize);

        nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
        decodability(nFold,:) = decodabilityLDA(nSessionData, trainingTargets, testTargets);
    end
    decodabilityAll(numRandPickUnits, :) = mean(decodability,1);
end
imagesc(nData_params.timeSeries, (1:maxRandPickUnits)*stepSize, decodabilityAll,[0.5 1]);
axis xy;
xlim([min(nData_params.timeSeries) max(nData_params.timeSeries)]);
ylim([stepSize maxRandPickUnits*stepSize])
gridxy ([nData_params.polein, nData_params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
box off;
hold off;
xlabel('Time (s)');
ylabel('# units');
colorbar
setPrint(8, 6, [PlotDir 'ModeledCollectedUnitsDecodability/CollectedUnitsDecodabilityDiffNumberUnitsFixedROCThres_' nData_name], 'pdf')


%%%%trialLDAEpoch
nDataSet = oldDataSet;
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
numFold             = 10;
selectedNeuronalIndex = nData_ActiveNeuronIndex';
selectedNeuronalIndex = selectedHighROCneurons(nDataSet, nData_params, ROCThres, selectedNeuronalIndex);
nDataSet              = nDataSet(selectedNeuronalIndex);
figure;
timePoints            = timePointTrialPeriod(nData_params.polein, nData_params.poleout, nData_params.timeSeries);
numPeriods            = length(timePoints) - 1;
decodability          = zeros(numFold, numPeriods, size(nDataSet(1).unit_yes_trial,2));
for nFold        = 1:numFold
    numUnits     = length(nDataSet);
    nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, numRandPickUnits)), totTargets, numTestTrials);
    nSessionData = permute(nSessionData,[1 3 2]);
    EpochIndex   = epochIndex(nData_params);
    EpochIndex   = EpochIndex(:,ones(1,numTrials))';
    decodability(nFold,:,:)    = decodabilitySliceDataTaperAccumulated(nSessionData, EpochIndex, 0, 0, numTestTrials, numPeriods);
end

subplot(10, 1, 1)
[~, maxIndex] = max(mean(decodability, 1), [], 2);
maxIndex = squeeze(maxIndex);
imagesc(1, nData_params.timeSeries, maxIndex', [1 4]);
axis off

subplot(10, 1, 2:9)
hold on
area(nData_params.timeSeries,squeeze(mean(decodability,1))','Edgecolor','none');
xlim([min(nData_params.timeSeries) max(nData_params.timeSeries)]);
ylim([0 1])
gridxy ([nData_params.polein, nData_params.poleout, 0],[0.5], 'Color','k','Linestyle','--','linewid', 0.5)
box off;
hold off;
ylabel('Decodability')
xlabel('Time (s)')

setPrint(8, 6, [PlotDir 'ModeledCollectedUnitsDecodabilityEpoch/CollectedUnitsDecodabilityEpoch_' nData_name], 'pdf')

%%%%trialLDAFixedNumberUnitsDiffROCThres
nDataSet = oldDataSet;
numFold             = 10;
numTrials           = 500;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
numRandPickUnits    = 100;
oldDataSet               = nDataSet;
figure;
hold on
for ROCThres          = 0.5:0.1:0.8;
    selectedNeuronalIndex = nData_ActiveNeuronIndex';
    selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, nData_params, ROCThres, selectedNeuronalIndex);
    nDataSet              = oldDataSet(selectedNeuronalIndex);
    decodability          = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));

    if numRandPickUnits< length(nDataSet)

        for nFold    = 1:numFold
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
            decodability(nFold,:) = decodabilityLDA(nSessionData, trainingTargets, testTargets);
        end

        plot(nData_params.timeSeries, mean(decodability,1));
    end
end
xlim([min(nData_params.timeSeries) max(nData_params.timeSeries)]);
ylim([0.5 1])
legend({'0.5', '0.6', '0.7', '0.8'},'location','southeast')
legend('boxoff')
gridxy ([nData_params.polein, nData_params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
box off;
hold off;
xlabel('Time (s)');
ylabel('Decodability');
setPrint(8, 6, [PlotDir 'ModeledCollectedUnitsDecodability/CollectedUnitsDecodabilityFixedNumberUnitsDiffROCThres_' nData_name], 'pdf')

%%%%trialLDALDAFixedNumberUnits
nDataSet = oldDataSet;
numRandPickUnits    = 100;
numTrials           = numRandPickUnits*10;
totTargets          = [true(numTrials/2,1); false(numTrials/2,1)];
ROCThres            = 0.5;
selectedNeuronalIndex = nData_ActiveNeuronIndex';
selectedNeuronalIndex = selectedHighROCneurons(nDataSet, nData_params, ROCThres, selectedNeuronalIndex);
nDataSet              = nDataSet(selectedNeuronalIndex);
numUnits              = length(nDataSet);
figure;
hold on
currRandPickUnits     = numRandPickUnits;
if currRandPickUnits>numUnits; currRandPickUnits = 100; end
nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, currRandPickUnits)), totTargets, numTrials);
nSessionData = normalizationDim(nSessionData, 2);
coeffs       = coeffLDA(nSessionData, totTargets);
imagesc(nData_params.timeSeries, nData_params.timeSeries, coeffs'*coeffs);
xlim([min(nData_params.timeSeries) max(nData_params.timeSeries)]);
ylim([min(nData_params.timeSeries) max(nData_params.timeSeries)]);
%     caxis([-1 1]);
colorbar
axis xy;
gridxy ([nData_params.polein, nData_params.poleout, 0],[nData_params.polein, nData_params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
box off;
hold off;
xlabel('Time (s)')
ylabel('Time (s)')
setPrint(8, 6, [PlotDir 'ModeledCollectedUnitsPCALDACorr/SimilarityLDALDA_' nData_name], 'pdf')

%%%%trialLDAPCAFixedNumberUnits
nDataSet = oldDataSet;
numRandPickUnits    = 400;
numTrials           = numRandPickUnits*3;
totTargets          = [true(numTrials/2,1); false(numTrials/2,1)];
ROCThres            = 0.5;
figure;
hold on
selectedNeuronalIndex = nData_ActiveNeuronIndex';
selectedNeuronalIndex = selectedHighROCneurons(nDataSet, nData_params, ROCThres, selectedNeuronalIndex);
nDataSet              = nDataSet(selectedNeuronalIndex);
numUnits              = length(nDataSet);
currRandPickUnits     = numRandPickUnits;
if currRandPickUnits>numUnits; currRandPickUnits = 100; end
nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, currRandPickUnits)), totTargets, numTrials);
nSessionData = normalizationDim(nSessionData, 2);
coeffPCAs    = coeffPCA(nSessionData);
coeffLDAs    = coeffLDA(nSessionData, totTargets);
imagesc(nData_params.timeSeries, nData_params.timeSeries, abs(coeffPCAs'*coeffLDAs));
xlim([min(nData_params.timeSeries) max(nData_params.timeSeries)]);
ylim([min(nData_params.timeSeries) max(nData_params.timeSeries)]);
%     caxis([0 0.6]);
colorbar
axis xy;
gridxy ([nData_params.polein, nData_params.poleout, 0],[nData_params.polein, nData_params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
box off;
hold off;
xlabel('Time (s)')
ylabel('Time (s)')
setPrint(8, 6, [PlotDir 'ModeledCollectedUnitsPCALDACorr/SimilarityLDAPCA_' nData_name], 'pdf')

%%%%trialLDAPCANoNormalizationFixedNumberUnits
nDataSet = oldDataSet;
numRandPickUnits    = 400;
numTrials           = numRandPickUnits*3;
totTargets          = [true(numTrials/2,1); false(numTrials/2,1)];
ROCThres            = 0.5;
figure;
hold on
selectedNeuronalIndex = nData_ActiveNeuronIndex';
selectedNeuronalIndex = selectedHighROCneurons(nDataSet, nData_params, ROCThres, selectedNeuronalIndex);
nDataSet              = nDataSet(selectedNeuronalIndex);
numUnits              = length(nDataSet);
currRandPickUnits     = numRandPickUnits;
if currRandPickUnits>numUnits; currRandPickUnits = 100; end
nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, currRandPickUnits)), totTargets, numTrials);
%     nSessionData = normalizationDim(nSessionData, 2);
coeffPCAs    = coeffPCA(nSessionData);
coeffLDAs    = coeffLDA(nSessionData, totTargets);
imagesc(nData_params.timeSeries, nData_params.timeSeries, abs(coeffPCAs'*coeffLDAs));
xlim([min(nData_params.timeSeries) max(nData_params.timeSeries)]);
ylim([min(nData_params.timeSeries) max(nData_params.timeSeries)]);
%     caxis([0 0.6]);
colorbar
axis xy;
gridxy ([nData_params.polein, nData_params.poleout, 0],[nData_params.polein, nData_params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
box off;
hold off;
xlabel('Time (s)')
ylabel('Time (s)')
setPrint(8, 6, [PlotDir 'CollectedUnitsPCALDACorr/SimilarityLDAPCANONormalization_' nData_name], 'pdf')

%%%%trialPCAAllUnits
nDataSet = oldDataSet;
numTrials      = 5000;
numComps       = 5;
trialType      = [true(numTrials, 1); false(numTrials, 1)];
numUnits           = length(nDataSet);
numT               = size(nDataSet(1).unit_yes_trial, 2);
pcaVar             = nan(numComps, numT);
evTrialType        = nan(numComps, numT);
firingRates        = generatePCAData(nDataSet, numTrials);


for nTime          = 1:numT    
    nFiringRates       = squeeze(firingRates(nTime, :, :))';    
    [~,score,latent,~,explained,~] = pca(nFiringRates, 'NumComponents', numComps);
    totVar            = sum(latent);
    pcaVar(:, nTime)  = explained(1:numComps)/sum(explained);
    evTrialType(:, nTime) = mean(score(1:numTrials, :)).^2/totVar;
end

save([TempDatDir 'PCATimeAllUnit_' nData_name '.mat'], 'pcaVar', 'evTrialType')

load([TempDatDir 'PCATimeAllUnit_' nData_name '.mat'], 'pcaVar', 'evTrialType')

figure;
hold on
imagesc(nData_params.timeSeries, 1:numComps, pcaVar(1:numComps, :))
axis xy;
xlim([min(nData_params.timeSeries) max(nData_params.timeSeries)]);
ylim([1 numComps])
gridxy ([nData_params.polein, nData_params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
box off;
hold off;
xlabel('Time (s)');
ylabel('Component index');
%     title('frac. PC Variance')
colorbar
setPrint(8, 6, [PlotDir 'ModeledCollectedUnitsPCA/CollectedUnitsPCAFracVar_' nData_name], 'pdf')

figure
hold on
imagesc(nData_params.timeSeries, 1:numComps, evTrialType(1:numComps, :))
axis xy;
xlim([min(nData_params.timeSeries) max(nData_params.timeSeries)]);
ylim([1 numComps])
gridxy ([nData_params.polein, nData_params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
box off;
hold off;
xlabel('Time (s)');
ylabel('Component index');
%     title('PC Trial Type EV')
colorbar

setPrint(8, 6, [PlotDir 'ModeledCollectedUnitsPCA/CollectedUnitsPCATrialEV_' nData_name], 'pdf')