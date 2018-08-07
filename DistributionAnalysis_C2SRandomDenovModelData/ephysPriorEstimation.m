%
% Compute ephys priors
%
% -------------------------------------------------------------------------
% version 1.0
%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

% all data is precomputed by code meanSynthEphysTraces.m

indexDatasets = [3, 4];
% 2: short Ca GP517
% 3: short Ca slow
% 4: short Ca slow virus
% 5: long Ca fast
% 6: long Ca slow

% load([TempDatDir DataSetList(1).name '.mat'])
% 
% numEphys  = length(nDataSet);
% numT      = length(DataSetList(1).params.timeSeries);
% actMat    = nan(length(nDataSet), numT * 2);
% 
% for nUnit = 1: length(nDataSet)
%     unitYesAct                = mean(nDataSet(nUnit).unit_yes_trial);
%     unitNoAct                 = mean(nDataSet(nUnit).unit_no_trial);
%     actMat(nUnit, 1:numT)     = unitYesAct;
%     actMat(nUnit, numT+1:end) = unitNoAct;
% end
% 
% ephysMat  = actMat;
% clear actMat nDataSet
% per_list           = 0.02:0.01:0.98;
% noise_factor_list  = sqrt(icdf('Exponential', per_list, 0.3527));
% params           = DataSetList(2).params;
% 
% sigma            = 0.15 / params.binsize; % 300 ms
% filterLength     = 11;
% filterStep       = linspace(-filterLength / 2, filterLength / 2, filterLength);
% filterInUse      = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
% filterInUse      = filterInUse / sum (filterInUse);

datDir        = '/Volumes/My Drive/ALM_Recording_Pole_Task_Svoboda_Lab/DatasetComparison/directDeconv/';

% for nData     = indexDatasets
%     for nTau  = 1:97
%         if exist([datDir DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'file')
%             load([datDir DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'],'spikeDataSet')
%             numSynEphys      = length(spikeDataSet);
%             ephysSimMat      = nan(numEphys, numSynEphys, length(noise_factor_list));
%             for nFF          = 1:length(noise_factor_list) 
%                 yesProfileMatrix = nan(length(spikeDataSet), numT);
%                 noProfileMatrix  = yesProfileMatrix;
%                 for nUnit        = 1:length(spikeDataSet)                    
%                     nUnitData        = spikeDataSet(nUnit).unit_yes_trial;
%                     yesUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
%                     nUnitData        = spikeDataSet(nUnit).unit_no_trial;
%                     noUnitData       = getGaussianPSTH (filterInUse, nUnitData, 2);
%                     mean_yesUnitData = mean(yesUnitData, 1);
%                     mean_noUnitData  = mean(noUnitData, 1);
%                     yesSize          = length(spikeDataSet(nUnit).unit_yes_trial_index);
%                     noSize           = length(spikeDataSet(nUnit).unit_no_trial_index);
%                     mean_rate        = min([mean_yesUnitData, mean_noUnitData]);
%                     baseline_rate    = 0;
%                     std_yesUnitData  = sqrt(mean_yesUnitData - mean_rate + baseline_rate);
%                     std_noUnitData   = sqrt(mean_noUnitData - mean_rate + baseline_rate);
%                     ffactor          = noise_factor_list(nFF);
%                     yesUnitData      = mean_yesUnitData + ones(size(yesUnitData,1), 1) * std_yesUnitData * ffactor .* randn(yesSize, numT);
%                     noUnitData       = mean_noUnitData + ones(size(noUnitData, 1), 1) * std_noUnitData * ffactor .* randn(noSize, numT);
%                     yesData          = mean(yesUnitData);
%                     noData           = mean(noUnitData);
%                     maxData          = max([yesData, noData]);
%                     minData          = min([yesData, noData]);
%                     rData            = (maxData - minData);
%                     yesData          = (yesData - minData)/(maxData - minData);
%                     noData           = (noData - minData)/(maxData - minData);
%                     yesProfileMatrix(nUnit, :)    = yesData;
%                     noProfileMatrix(nUnit, :)     = noData;
%                 end
%                 actMat                 = [yesProfileMatrix, noProfileMatrix];
%                 ephysSimMat(:, :, nFF) = corr(ephysMat', actMat', 'type', 'Spearman');
%             end  
%             save([datDir DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '_ephysSimMat.mat'], 'ephysSimMat');
%             clear spikeDataSet ephysSimMat
%         end
%     end    
% end

for nData           = indexDatasets
    ephysSimMatSum  = 0;
    for nTau        = 1:97
        if exist([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'file')
            load([datDir DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '_ephysSimMat.mat'], 'ephysSimMat')
            ephysSimMat(ephysSimMat<0) = 0;
            ephysSimMatSum             = ephysSimMatSum + squeeze(sum(ephysSimMat, 3));
            clear ephysSimMat
        end
    end   
    figure;
    imagesc(ephysSimMatSum)
    save([TempDatDir 'directDeconv/' DataSetList(nData).name '_ephysSimSumMat.mat'], 'ephysSimMatSum');
end