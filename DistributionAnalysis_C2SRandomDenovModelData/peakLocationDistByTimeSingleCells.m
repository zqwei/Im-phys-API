%
% Compute distribution of selectivity for each single neuron
% 
% Here we drop the test of parameter of tau_r
% 
% 
% -------------------------------------------------------------------------
% version 1.0
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

% all data is precomputed by code meanSynthEphysTraces.m

indexDatasets = [2, 3, 4, 5, 6];
% 2: short Ca GP517
% 3: short Ca slow
% 4: short Ca slow virus
% 5: long Ca fast
% 6: long Ca slow

% params           = DataSetList(2).params;
% 
% sigma            = 0.15 / params.binsize; % 300 ms
% filterLength     = 11;
% filterStep       = linspace(-filterLength / 2, filterLength / 2, filterLength);
% filterInUse      = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
% filterInUse      = filterInUse / sum (filterInUse); 
% 
% per_list           = 0.02:0.01:0.98;
% noise_factor_list  = sqrt(icdf('Exponential', per_list, 0.3527));
% timeTag            = 8:60;
% 
% for nData     = indexDatasets
%     for nTau  = 1:99
%         if exist([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'file')
%             load([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'],'spikeDataSet')
%             peakiness        = nan(99, 1);
%             for nFF          = 1:length(noise_factor_list) 
%                 numTimeBin       = size(spikeDataSet(1).unit_yes_trial, 2);
%                 yesProfileMatrix = nan(length(spikeDataSet), numTimeBin);
%                 noProfileMatrix  = yesProfileMatrix;
%                 positivePeak     = false(length(spikeDataSet));
%                 for nUnit        = 1:length(spikeDataSet)
%                 nUnitData        = spikeDataSet(nUnit).unit_yes_trial;
%                 yesUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
%                 nUnitData        = spikeDataSet(nUnit).unit_no_trial;
%                 noUnitData       = getGaussianPSTH (filterInUse, nUnitData, 2);
%                 mean_yesUnitData = mean(yesUnitData, 1);
%                 mean_noUnitData  = mean(noUnitData, 1);
%                 yesSize          = length(spikeDataSet(nUnit).unit_yes_trial_index);
%                 noSize           = length(spikeDataSet(nUnit).unit_no_trial_index);
%                 numT             = size(nUnitData, 2);
%                 mean_rate        = min([mean_yesUnitData, mean_noUnitData]);
%                 baseline_rate    = 0;
%                 std_yesUnitData  = sqrt(mean_yesUnitData - mean_rate + baseline_rate);
%                 std_noUnitData   = sqrt(mean_noUnitData - mean_rate + baseline_rate);
%                 ffactor          = noise_factor_list(nFF);
%                 yesUnitData      = mean_yesUnitData + ones(size(yesUnitData,1), 1) * std_yesUnitData * ffactor .* randn(yesSize, numT);
%                 noUnitData       = mean_noUnitData + ones(size(noUnitData, 1), 1) * std_noUnitData * ffactor .* randn(noSize, numT);
%                 yesData          = mean(yesUnitData);
%                 noData           = mean(noUnitData);
%                 maxData          = max([yesData, noData]);
%                 minData          = min([yesData, noData]);
%                 rData            = (maxData - minData);
%                 yesData          = (yesData - minData)/(maxData - minData);
%                 noData           = (noData - minData)/(maxData - minData);
%                 yesProfileMatrix(nUnit, :)    = yesData;
%                 noProfileMatrix(nUnit, :)     = noData;
%                 positivePeak(nUnit)       = mean(yesData(1:8)) <= mean(yesData(9:47)) ...
%                                            || mean(noData(1:8)) <= mean(noData(9:47));
%                 end
%                 actMat        = [yesProfileMatrix, noProfileMatrix];
%                 actMat        = actMat(positivePeak, :);
%                 [~, maxId]    = max(actMat, [], 2);
%                 countMaxId    = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
%                 peakiness(nFF) = std(countMaxId([timeTag, timeTag+77]));
%             end  
%             save([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '_peak.mat'], 'peakiness');
%             clear spikeDataSet peakiness
%         end
%     end    
% end

per_list           = 0.02:0.01:0.98;

for nData     = [3, 4] %indexDatasets
    s_mat     = nan(97, 97);
    for nTau  = 1:97
        if exist([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'file')
            load([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '_peak.mat'], 'peakiness')
            s_mat(nTau, :) = peakiness(1:97);
        end
    end    
    disp(mean(s_mat(:)))
    disp(std(s_mat(:)))
    figure;
    imagesc(per_list*100, per_list*100, s_mat-1.27, [-0.9 -0.7]);
    colormap(gray)
    colorbar
    xlabel('percentile noise level')
    ylabel('percentile \tau_d')
    axis xy
    colorbar
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'DistributionAnanlysis/Peakiness_' DataSetList(nData).name], 'pdf')
end