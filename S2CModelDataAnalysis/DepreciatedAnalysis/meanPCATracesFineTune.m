%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load('ParamsFitCells_S2CModel_Sim.mat');
S2Cparams   = params;
clear params;
idTauD      = 1;
idN         = 1;
idK         = 1;
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

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numComps       = 3;

load ([TempDatDir 'DataListShuffle.mat']);
nData = 1;
load([TempDatDir DataSetList(nData).name '.mat']);
% logPValueEpoch     = getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
% unitGroup          = plotTtestLogPSpikeEpoch (logPValueEpoch);
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
params.Fm              = ones(length(nDataSet), 1)*21.3240;
params.K               = ones(length(nDataSet), 1)*S2Cparams(1).K;
params.n               = ones(length(nDataSet), 1)*S2Cparams(1).n;
% params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  S2Cparams(12).tau_r + S2Cparams(1).tau_r;
% params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  S2Cparams(12).tau_d + S2Cparams(idTauD).tau_d;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  0 + S2Cparams(1).tau_r;
params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  0 + S2Cparams(idTauD).tau_d;
params.Fm              = params.Fm;
params.tau_r           = params.tau_r;
params.tau_d           = params.tau_d*1.0;
params.intNoise        = 1.5;
params.extNoise        = 0;
% nDataSet               = getFakeCaImagingData(nDataSet(unitGroup>0), params);
oldDataSet             = getFakeCaImagingData(nDataSet, params);
% firingRates        = generateDPCADataLinearData(nDataSet, numTrials);
% firingRatesAverage = nanmean(firingRates, ndims(firingRates));
% firingRatesAverage = [squeeze(firingRatesAverage(:, 1, :)), squeeze(firingRatesAverage(:, 2, :));];
% [~,score,~]        = pca(firingRatesAverage', 'NumComponents', numComps);
% scoreLinear        = score(:, 1);
% % scoreLinear        = score(1:77, 1);
% scoreLinear        = scoreLinear - min(scoreLinear);
% % figure('Visible', 'on')
% % plot(scoreLinear)
% k  = 9;
% plotKMeansLinearFilter(nDataSet, k, params);
% plotKMeans(nDataSet, k, params);

% oldDataSet         = nDataSet;
nDataSet           = generateNonLinearDataFromLinear(oldDataSet, params);
firingRates        = generateDPCAData(nDataSet(contraIndex), numTrials);
firingRatesAverage = nanmean(firingRates, ndims(firingRates));
% pcaX               = firingRatesAverage(:,:);
% firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
% firingRatesAverage = bsxfun(@rdivide, firingRatesAverage, std(pcaX,[],2));
firingRatesAverage = [squeeze(firingRatesAverage(:, 1, :)), squeeze(firingRatesAverage(:, 2, :));];
[coeffs,score,~]        = pca(firingRatesAverage', 'NumComponents', numComps);
scoreS2CLinear     = score(:, 1);
% scoreS2CLinear        = score(1:77, 1);
scoreS2CLinear     = scoreS2CLinear - min(scoreS2CLinear);
% plotKMeans(nDataSet, k, params);
figure('Visible', 'on')
for ii = 1:numComps
    subplot(1, numComps, ii)
    plot([score(1:77, ii),score(78:end, ii)])
end
figure('Visible', 'on')
plot([tanh(mean(firingRatesAverage(:, 1:77))')+0.5,tanh(mean(firingRatesAverage(:,78:end))')+0.5])

% figure;
% bar(coeffs(:,1))
% 
% figure;
% hist(coeffs(:,1), 100)
% 
% figure;
% plot([firingRatesAverage(398, 1:77)',firingRatesAverage(398,78:end)'])
% 



% load ([TempDatDir 'DataListShuffle.mat']);
% nData = 4;
% load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat']);
% firingRates        = generateDPCAData(nDataSet, numTrials);
% firingRatesAverage = nanmean(firingRates, ndims(firingRates));
% pcaFiringRatesAverage = zeros(numComps, 2, 77);
% firingRatesAverage = [squeeze(firingRatesAverage(:, 1, :)), squeeze(firingRatesAverage(:, 2, :));];
% [~,score,~]        = pca(firingRatesAverage', 'NumComponents', numComps);
% scoreNonLinear     = -score(:, 1);
% % scoreNonLinear     = -score(1:77, 1);
% scoreNonLinear     = scoreNonLinear - min(scoreNonLinear);
% % valid_index        = scoreNonLinear>0;
% 
% hill_fit = @(b,x)  b(1).*x.^b(3)./(b(2)+x.^b(3));
% b0 = [15; 20; 2];
% B = lsqcurvefit(hill_fit, b0, scoreLinear(1:77), scoreNonLinear(1:77));
% 
% figure('Visible', 'on')
% hold on
% % plot(scoreLinear(1:77), scoreNonLinear(1:77), 'ok')
% % plot(scoreLinear(78:end), scoreNonLinear(78:end), 'or')
% % linearVec          = linspace(0, 50);
% % plot(linearVec, hill_fit(B,linearVec), '-r')
% plot([scoreS2CLinear, scoreNonLinear])