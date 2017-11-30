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



% % % computing cell type and save
% per_list               = 0.02:0.01:0.98;
% noise_factor_list      = sqrt(icdf('Exponential', per_list, 0.3527));
% 
% for nData     = indexDatasets
%     params    = DataSetList(nData).params;
%     for nTau  = 1:99
%         if exist([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'file')
%             load([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'],'spikeDataSet')
%             cellTypeMat       = nan(size(spikeDataSet, 1), 99);
%             for nFF           = 1:length(noise_factor_list) 
%               unitGroup       = getLogPValueTscoreSpikeTimeAve(spikeDataSet, params, noise_factor_list(nFF));
%               cellTypeMat(:, nFF) = unitGroup;
%             end  
%             save([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '_cell_type.mat'], 'cellTypeMat');
%             clear spikeDataSet cellTypeMat
%         end
%     end    
% end


% % % plot saved file without ephys prior data
% per_list               = 0.02:0.01:0.98;
% 
% for nData     = indexDatasets
%     monoCell  = nan(97, 97);
%     multiCell = nan(97, 97);
%     validTau  = false(97, 1);
%     for nTau  = 1:99
%         if exist([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'file')
%             load([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '_cell_type.mat']);
%             monoCell(nTau, :)  = mean(cellTypeMat(:, 1:97) == 1);
%             multiCell(nTau, :) = mean(cellTypeMat(:, 1:97) == 2);
%             validTau(nTau)     = true;
%         end
%     end    
%     
%     figure;
%     subplot(1, 2, 1)
%     imagesc(per_list(validTau)*100, per_list*100, monoCell(validTau, :), [0 1])
%     axis xy
%     xlabel('Noise level (%)')
%     ylabel('\tau_d (%)')
%     title('Frac Monophasic neuron')
%     subplot(1, 2, 2)
%     imagesc(per_list(validTau), per_list, multiCell(validTau, :), [0 1])
%     axis xy
%     xlabel('Noise level (%)')
%     ylabel('\tau_d (%)')
%     title('Frac Multiphasic neuron')
%     setPrint(8*2, 6, [PlotDir 'DistributionAnanlysis/Cell_type_' DataSetList(nData).name])
% end


% % % plot saved file with ephys prior data


per_list      = 0.02:0.01:0.98;

for nData     = indexDatasets
    monoCell  = nan(97, 97);
    multiCell = nan(97, 97);
    validTau  = false(97, 1);
    for nTau  = 1:99
        if exist([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'file')
            load([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '_cell_type.mat']);
            monoCell(nTau, :)  = mean(cellTypeMat(:, 1:97) == 1);
            multiCell(nTau, :) = mean(cellTypeMat(:, 1:97) == 2);
            validTau(nTau)     = true;
        end
    end    
    
    figure;
    subplot(1, 2, 1)
    imagesc(per_list(validTau)*100, per_list*100, monoCell(validTau, :), [0 1])
    axis xy
    xlabel('Noise level (%)')
    ylabel('\tau_d (%)')
    title('Frac Monophasic neuron')
    subplot(1, 2, 2)
    imagesc(per_list(validTau), per_list, multiCell(validTau, :), [0 1])
    axis xy
    xlabel('Noise level (%)')
    ylabel('\tau_d (%)')
    title('Frac Multiphasic neuron')
    setPrint(8*2, 6, [PlotDir 'DistributionAnanlysis/Cell_type_' DataSetList(nData).name])
end