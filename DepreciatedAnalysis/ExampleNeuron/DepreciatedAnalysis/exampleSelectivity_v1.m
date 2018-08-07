addpath('../Func');
setDir;

nCell = 290;
load ([TempDatDir 'DataListShuffle.mat']);
nData = 1; % plot raster and psth
load([TempDatDir DataSetList(nData).name '.mat'])
spikeDataSet = nDataSet;   
params = DataSetList(nData).params;
load ([TempDatDir 'DataListS2CModel.mat']);    
nData = 4; % plot s2c model
load([TempDatDir DataSetList(nData).name '.mat'])
s2cDataSet   = nDataSet;

figure;
% spike
subplot(1,2, 1)
sigma                         = 0.15 / params.binsize; % 200 ms
filterLength                  = 11;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse); 

nUnitData        = spikeDataSet(nCell).unit_yes_trial;
yesUnitData      = nUnitData;%getGaussianPSTH (filterInUse, nUnitData, 2);
nUnitData        = spikeDataSet(nCell).unit_no_trial;
noUnitData       = nUnitData;%getGaussianPSTH (filterInUse, nUnitData, 2);
% [~, p, ~, stats]           = ttest2(yesUnitData, noUnitData);
% p(p<exp(-10))    = exp(-10);
d_prime          = (mean(yesUnitData) - mean(noUnitData))./sqrt((var(yesUnitData)+var(noUnitData))/2);
d_prime          = getGaussianPSTH (filterInUse, d_prime, 2);
hold on
plot(params.timeSeries, d_prime, 'Color','k','Linestyle','-','linewid', 1.0)
ylim([0 3])
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5);
hold off;
ylabel('two-sample t-score');
xlabel('Time (s)');
xlim([params.timeSeries(1) params.timeSeries(end)]);
ylim([0 3])
set(gca, 'TickDir', 'out')
%
subplot(1, 2, 2)
nUnitData        = s2cDataSet(nCell).unit_yes_trial;
yesUnitData      = nUnitData;%getGaussianPSTH (filterInUse, nUnitData, 2);
nUnitData        = s2cDataSet(nCell).unit_no_trial;
noUnitData       = nUnitData;%getGaussianPSTH (filterInUse, nUnitData, 2);
% [~, p, ~, stats]           = ttest2(yesUnitData, noUnitData);
% p(p<exp(-10))    = exp(-10);
d_prime          = (mean(yesUnitData) - mean(noUnitData))./sqrt((var(yesUnitData)+var(noUnitData))/2);
d_prime          = getGaussianPSTH (filterInUse, d_prime, 2);
hold on
plot(params.timeSeries, d_prime, 'Color','k','Linestyle','-','linewid', 1.0)
ylim([0 3])
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5);
hold off;
ylabel('two-sample t-score');
xlabel('Time (s)');
xlim([params.timeSeries(1) params.timeSeries(end)]);
ylim([0 3])
set(gca, 'TickDir', 'out')
setPrint(8*2, 6*1, ['SingleUnitsTscoreExampleNeuron'])

% function plotRaster(spikeDataSet, params)
%     spkTimes{1}    = spikeDataSet.unit_yes_trial_spk_time;
%     spkTimes{2}    = spikeDataSet.unit_no_trial_spk_time;
%     hold on;
%     color_index    = [0.7  0 0; 0 0 0.7];
%     spkLoc         = 0;
%     for nPlot            = 1:2
%         hold on;
%         for nTrial       = 1:20
%             spkLoc       = spkLoc + 1;
%             spikeTimes   = spkTimes{nPlot}{nTrial};
%             spikeTrial   = ones(length(spikeTimes), 1) * spkLoc;
%             plot(spikeTimes, spikeTrial, '.', 'color', color_index(nPlot, :));
%         end
%         spkLoc     = spkLoc + 1;%length(spkTimes{nPlot}) + 3;
%     end
% 
%     gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0);
%     hold off;
%     ylim([1 41])
%     xlim([params.timeSeries(1) params.timeSeries(end)]);
%     axis off
% end
% 
% function plotDff(spikeDataSet, params)
%     sigma                         = 0.15 / params.binsize; % 200 ms
%     filterLength                  = 11;
%     filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
%     filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
%     filterInUse                   = filterInUse / sum (filterInUse); 
% 
%     nUnitData        = spikeDataSet.unit_yes_trial;
%     yesUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
%     nUnitData        = spikeDataSet.unit_no_trial;
%     noUnitData       = getGaussianPSTH (filterInUse, nUnitData, 2);
%     cmin             = min([mean(yesUnitData), mean(noUnitData)]);
%     cmax             = max([mean(yesUnitData), mean(noUnitData)]);
%     
% 
% 
%     hold on
%     actMat = nan(41, size(spikeDataSet.unit_yes_trial, 2));
%     actMat(1:20, :) = yesUnitData(1:20, :);
%     actMat(22:end, :) = noUnitData(1:20, :);
%     imagesc(params.timeSeries, 1:41, actMat);
%     gridxy ([params.polein, params.poleout, 0],[], 'Color','w','Linestyle','--','linewid', 1.0);
%     hold off;
%     ylim([1 41])
%     colormap(gray)
%     
%     caxis([cmin, cmax]);
% %     colorbar
%     xlim([params.timeSeries(1) params.timeSeries(end)]);
%     axis off
% end
% 
% 
% function plotPSTH(spikeDataSet, params, ylabelName)
%     sigma                         = 0.15 / params.binsize; % 200 ms
%     filterLength                  = 11;
%     filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
%     filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
%     filterInUse                   = filterInUse / sum (filterInUse); 
% 
%     color_index    = [0.7  0 0; 0 0 0.7];
%     hold on;
%     nUnitData        = spikeDataSet.unit_yes_trial;
%     nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
%     shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
%         std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
%         {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
%     nUnitData        = spikeDataSet.unit_no_trial;
%     nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
%     shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
%         std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
%         {'k-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
%     gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
%     hold off;
%     ylabel(ylabelName);
%     xlabel('Time (s)');
%     xlim([params.timeSeries(1) params.timeSeries(end)]);
%     set(gca, 'TickDir', 'out')
% end
% 
% 
% function meanDataSet = meanData(nDataSet)
%     meanDataSet = nan(length(nDataSet), size(nDataSet(1).unit_yes_trial, 2)*2);
%     for nCell = 1:length(nDataSet)
%         meanDataSet(nCell, :) = [mean(nDataSet(nCell).unit_yes_trial), mean(nDataSet(nCell).unit_no_trial)];
%     end
% end
% 
% function mIndex = findSimilarCellToS2CModel(nCellS2CDataSet, caDataSet)
%     corrDat     = corr(nCellS2CDataSet', caDataSet');
%     [~, mIndex] = max(corrDat);
%     mIndex      = mIndex(1);
% end