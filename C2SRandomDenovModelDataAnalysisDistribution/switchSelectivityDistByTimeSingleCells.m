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

indexDatasets = [2, 3, 5, 6];
% 2: short Ca GP517
% 3: short Ca slow
% 5: long Ca fast
% 6: long Ca slow




per_list               = 0.02:0.01:0.98;
noise_factor_list      = sqrt(icdf('Exponential', per_list, 0.3527));

for nData     = indexDatasets
    for nTau  = 1:99
        if exist([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'file')
            load([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'],'spikeDataSet')
            cellTypeMat       = nan(size(spikeDataSet, 1), 99);
            for nFF           = 1:length(noise_factor_list) 
              unitGroup       = getLogPValueTscoreSpikeTimeAve(spikeDataSet, params, noise_factor_list(nFF));
              cellTypeMat(:, nFF) = unitGroup;
            end  
            save([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '_cell_type.mat'], 'cellTypeMat');
        end
    end    
end