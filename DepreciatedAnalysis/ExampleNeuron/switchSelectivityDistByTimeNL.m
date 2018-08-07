%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

cmap = [0.8000    0.8000    0.8000;
       1.0000    0.6000         0;
       0    0.8000         0];

nData   = 1;
load([TempDatDir DataSetList(nData).name '.mat'])
spikeDataSet = nDataSet;   


% timeLength       = size(spikeDataSet(1).unit_yes_trial, 2);
fanoFactorMat    = nan(length(nDataSet),  2);

for nUnit        = 1:length(nDataSet)
    yesfanoFactor           = var(mean(spikeDataSet(nUnit).unit_yes_trial, 2))./mean(mean(spikeDataSet(nUnit).unit_yes_trial, 2));
    nofanoFactor            = var(mean(spikeDataSet(nUnit).unit_no_trial, 2))./mean(mean(spikeDataSet(nUnit).unit_no_trial, 2));
    fanoFactorMat(nUnit, :) = [yesfanoFactor, nofanoFactor];
end

figure
hist(fanoFactorMat(:), 100,'facecolor','k')
[f, xi] = ksdensity(fanoFactorMat(:));
hold on
plot(xi, f/max(f)*144, 'r')
xlim([0 8])
setPrint(8, 6, 'Fanofactor')

mean_fr            = 0.0;

load('DeconvGPResults/S2CC2SNLRandomLinearAve.mat', 'nDataSet')
fracNeuron = [];
for ffactor = 1.0:0.05:4.4
    unitGroup = getLogPValueTscoreSpikeTimeAve(nDataSet, DataSetList(nData).params, ffactor);    
    sizeGroup = histcounts(unitGroup, 0:3);
    fracNeuron = [fracNeuron; sizeGroup(2:3)/sum(sizeGroup)];
end

figure;
plot(1.0:0.05:4.4, fracNeuron, '-o')
xlabel('variability -- Fano factor')
ylabel('Frac. neuron')
legend('Mono.', 'Multi.')
xlim([0.99 4.41])
ylim([0 1])
setPrint(8, 6, 'Fano factor S2C')

% fanoFactorConstant = 0.8;
% 
% for nUnit        = 1:length(nDataSet)
%     mean_rate    = min([nDataSet(nUnit).decv_yes_trial(:); nDataSet(nUnit).decv_no_trial(:)]);
%     if size(nDataSet(nUnit).decv_yes_trial, 1) > 1
%         nDataSet(nUnit).unit_yes_trial      = nDataSet(nUnit).decv_yes_trial;
% %         nDataSet(nUnit).unit_yes_trial(nDataSet(nUnit).unit_yes_trial<0) = 0;
%     else
%         fanoFactor                          = std(spikeDataSet(nUnit).unit_yes_trial)./mean(spikeDataSet(nUnit).unit_yes_trial);
%         stdTrial                            = sqrt(nDataSet(nUnit).decv_yes_trial - mean_rate + mean_fr) * fanoFactorConstant;% .* fanoFactor;
%         yesNoise                            = nan(size(spikeDataSet(nUnit).unit_yes_trial));
%         for nTrial = 1:size(spikeDataSet(nUnit).unit_yes_trial, 1)
%             yesNoise(nTrial, :)             = randn(1, size(spikeDataSet(nUnit).unit_yes_trial, 2)) .* stdTrial;
%         end
%         yesNoise(isnan(yesNoise))           = 0;
%         nDataSet(nUnit).unit_yes_trial      = nDataSet(nUnit).decv_yes_trial + yesNoise;
% %         nDataSet(nUnit).unit_yes_trial(nDataSet(nUnit).unit_yes_trial<0) = 0;
%     end
%         
%     if size(nDataSet(nUnit).decv_no_trial, 1) > 1
%         nDataSet(nUnit).unit_no_trial       = nDataSet(nUnit).decv_no_trial;
% %         nDataSet(nUnit).unit_no_trial(nDataSet(nUnit).unit_no_trial<0) = 0;
%     else
%         fanoFactor                          = std(spikeDataSet(nUnit).unit_no_trial)./mean(spikeDataSet(nUnit).unit_no_trial);
%         stdTrial                            = sqrt(nDataSet(nUnit).decv_no_trial - mean_rate + mean_fr) * fanoFactorConstant;% .* fanoFactor;
%         noNoise                             = nan(size(spikeDataSet(nUnit).unit_no_trial));
%         for nTrial = 1:size(spikeDataSet(nUnit).unit_no_trial, 1)
%             noNoise(nTrial, :)              = randn(1, size(spikeDataSet(nUnit).unit_no_trial, 2)) .* stdTrial;
%         end
%         noNoise(isnan(noNoise))             = 0;
%         nDataSet(nUnit).unit_no_trial       = nDataSet(nUnit).decv_no_trial + noNoise;
% %         nDataSet(nUnit).unit_no_trial(nDataSet(nUnit).unit_no_trial<0) = 0;
%     end
% end

load('DeconvResults/S2CC2SNLRandomNLSingleTrial.mat', 'nDataSet')
unitGroup = getLogPValueTscoreSpikeTimeAve(nDataSet, DataSetList(nData).params, 2.5); 
sizeGroup = histcounts(unitGroup, 0:3);

figure;
groupNames      = {'Non.', 'Homo.', 'Dynamical'};
pie(sizeGroup)
colormap(cmap)
set(gca, 'TickDir', 'out')
setPrint(8, 6, 'SingleUnitsTscoreTime')


addpath('../Func');
setDir;
load([TempDatDir 'Shuffle_Spikes.mat'])
unitSpikeGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
spikeDataSet    = nDataSet;
sigma                         = 0.15 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 10;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);
% load ([TempDatDir 'DataListS2CModel.mat']);

% nData   = 1;

% load('S2CC2S.mat', 'nDataSet')
%     
% for nUnit        = 1:length(nDataSet)
%     if size(nDataSet(nUnit).decv_yes_trial, 1) > 1
%         nDataSet(nUnit).unit_yes_trial      = nDataSet(nUnit).decv_yes_trial;
% %         nDataSet(nUnit).unit_yes_trial(nDataSet(nUnit).unit_yes_trial<0) = 0;
%     else
%         fanoFactor                          = std(spikeDataSet(nUnit).unit_yes_trial)./mean(spikeDataSet(nUnit).unit_yes_trial);
%         stdTrial                            = sqrt(nDataSet(nUnit).decv_yes_trial - mean_rate + mean_fr) * fanoFactorConstant;% .* fanoFactor;
%         yesNoise                            = nan(size(spikeDataSet(nUnit).unit_yes_trial));
%         for nTrial = 1:size(spikeDataSet(nUnit).unit_yes_trial, 1)
%             yesNoise(nTrial, :)             = randn(1, size(spikeDataSet(nUnit).unit_yes_trial, 2)) .* stdTrial;
%         end
%         yesNoise(isnan(yesNoise))           = 0;
%         nDataSet(nUnit).unit_yes_trial      = nDataSet(nUnit).decv_yes_trial + yesNoise;
% %         nDataSet(nUnit).unit_yes_trial(nDataSet(nUnit).unit_yes_trial<0) = 0;
%     end
%         
%     if size(nDataSet(nUnit).decv_no_trial, 1) > 1
%         nDataSet(nUnit).unit_no_trial       = nDataSet(nUnit).decv_no_trial;
% %         nDataSet(nUnit).unit_no_trial(nDataSet(nUnit).unit_no_trial<0) = 0;
%     else
%         fanoFactor                          = std(spikeDataSet(nUnit).unit_no_trial)./mean(spikeDataSet(nUnit).unit_no_trial);
%         stdTrial                            = sqrt(nDataSet(nUnit).decv_no_trial - mean_rate + mean_fr) * fanoFactorConstant;% .* fanoFactor;
%         noNoise                             = nan(size(spikeDataSet(nUnit).unit_no_trial));
%         for nTrial = 1:size(spikeDataSet(nUnit).unit_no_trial, 1)
%             noNoise(nTrial, :)              = randn(1, size(spikeDataSet(nUnit).unit_no_trial, 2)) .* stdTrial;
%         end
%         noNoise(isnan(noNoise))           = 0;
%         nDataSet(nUnit).unit_no_trial       = nDataSet(nUnit).decv_no_trial + noNoise;
% %         nDataSet(nUnit).unit_no_trial(nDataSet(nUnit).unit_no_trial<0) = 0;
%     end
% end
% unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params); 

load('DeconvGPResults/S2CC2SNLRandomLinearAve.mat', 'nDataSet')
unitGroup = getLogPValueTscoreSpikeTimeAve(nDataSet, DataSetList(nData).params, 2.5); 

load('DeconvGPResults/S2CC2SNLPreciseSingleTrial.mat', 'nDataSet')
spkDataSet = nDataSet;
for nCell  = 1:length(nDataSet)
    spkDataSet(nCell).unit_yes_trial = spkDataSet(nCell).decv_yes_trial;
    spkDataSet(nCell).unit_no_trial  = spkDataSet(nCell).decv_no_trial;
end
unitGroup = getLogPValueTscoreSpikeTime(spkDataSet, DataSetList(nData).params); 

load('DeconvGPResults/S2CC2SMCMCSingleTrial_FineTunedModeled_GP43.mat')
spkDataSet = nDataSet;
for nCell  = 1:length(nDataSet)
    spkDataSet(nCell).unit_yes_trial = spkDataSet(nCell).mcmc_yes_trial;
    spkDataSet(nCell).unit_no_trial  = spkDataSet(nCell).mcmc_no_trial;
end
unitGroup = getLogPValueTscoreSpikeTime(spkDataSet, DataSetList(nData).params); 

    
groupCounts = zeros(3, 3);

for nGroup  = 1:3
    for mGroup = 1:3
        groupCounts(nGroup, mGroup) = sum(unitGroup == mGroup-1 & unitSpikeGroup == nGroup-1);
    end
end

groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));

figure;
for nPlot = 1:3
    subplot(1, 3, nPlot)
    hold on
    bar(1:3, groupPerCounts(nPlot, :), 'edgecolor', 'none');
    ylim([0 1])
    box off
    xlabel('cell type')
    ylabel('frac. cell')
    colormap(cmap)
    xlim([0.5 3.5])
    set(gca, 'TickDir', 'out')
end
setPrint(8*3, 6, ['SingleUnitsTscoreTimeTrans_S2CC2S_' DataSetList(nData).name])

% groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 1));
% plot(diag(groupPerCounts), '-o')
% ylim([0 1])
% xlim([0.5 3.5])
% xlabel('cell type')
% ylabel('True positive rate')
% setPrint(8, 6, ['SingleUnitsTscoreTimeTrans_S2CC2STPR_' DataSetList(nData).name])



cmap = [0.8000    0.8000    0.8000; 0.4000    0.4000    0.4000;
        1.0000    0.6000         0; 1.0000    0.2000         0;
             0    0.8000         0;      0    0.4000         0];

truePos = diag(groupCounts)';
sizeGroup = sum(groupCounts, 1);
trueNeg = sizeGroup - truePos;
sizeGroup = [truePos; trueNeg];

figure;
groupNames      = {'Non.', 'Homo.', 'Dynamical'};
pie(sizeGroup(:))
colormap(cmap)
set(gca, 'TickDir', 'out')
setPrint(8, 6, 'SingleUnitsTscoreTimeTruePos')



% params = DataSetList(nData).params;
% nCell = 511;% 443; %511; % 511 % 662
% figure
% subplot(1, 2, 1)
% hold on
% nUnitData        = spikeDataSet(nCell).unit_yes_trial;
% nUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
% shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
%     std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
%     {'b-', 'linewid', 1.0}, 0.5);
% nUnitData        = spikeDataSet(nCell).unit_no_trial;
% nUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
% shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
%     std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
%     {'r-', 'linewid', 1.0}, 0.5);
% 
% subplot(1, 2, 2)
% hold on
% nUnitData        = nDataSet(nCell).unit_yes_trial;
% nUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
% shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
%     std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
%     {'b-', 'linewid', 1.0}, 0.5);
% nUnitData        = nDataSet(nCell).unit_no_trial;
% nUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
% shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
%     std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
%     {'r-', 'linewid', 1.0}, 0.5);

% close all
