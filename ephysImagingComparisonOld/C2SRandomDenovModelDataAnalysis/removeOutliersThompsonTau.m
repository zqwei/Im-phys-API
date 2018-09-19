%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
caDataSetList = DataSetList;
load ([TempDatDir 'DataListC2SRandomDeconvModel.mat']);
spikeDataSetList = DataSetList;

thres          = 0.3;

for nData              = [1 2]
    load([TempDatDir caDataSetList(nData + 2).name '_withOLRemoval.mat'])
%     caDataSet          = nDataSet;
%     params             = caDataSetList(nData).params;
%     timePeriod         = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
%     neuronRemoveList   = false(length(nDataSet), 1);
    load([TempDatDir spikeDataSetList(nData).name '.mat'])
    spikeDataSet       = nDataSet;
        
%     for nUnit          = 1:length(caDataSet)
%         caCellData     = caDataSet(nUnit);
%         spikeCellData  = spikeDataSet(nUnit);
%         idx_yes_remove = mean(isnan(caCellData.unit_yes_trial_removeoutlier), 2)>thres;
%         idx_no_remove  = mean(isnan(caCellData.unit_no_trial_removeoutlier), 2)>thres; 
%         if sum(idx_yes_remove)>0 || sum(idx_no_remove)>0 
%             if sum(~idx_yes_remove)>20 && sum(~idx_no_remove)>20
%                 spikeCellData.unit_yes_trial       = spikeCellData.unit_yes_trial(~idx_yes_remove, :);
%                 spikeCellData.unit_yes_trial_index = spikeCellData.unit_yes_trial_index(~idx_yes_remove);
%                 spikeCellData.unit_no_trial        = spikeCellData.unit_no_trial(~idx_no_remove, :);
%                 spikeCellData.unit_no_trial_index  = spikeCellData.unit_no_trial_index(~idx_no_remove);
%                 
%                 if ~ttest2(mean(caCellData.unit_yes_trial(:, timePeriod(1):timePeriod(2))), mean(caCellData.unit_no_trial(:, timePeriod(1):timePeriod(2))))
%                     spikeDataSet(nUnit) = spikeCellData;
%                 else
%                     neuronRemoveList(nUnit) = true;
%                 end
%             else
%                 neuronRemoveList(nUnit) = true;
%             end
%         end  
%     end
    spikeDataSet = spikeDataSet(~neuronRemoveList);
    nDataSet     = spikeDataSet;
    save([TempDatDir spikeDataSetList(nData).name '_withOLRemoval.mat'], 'nDataSet', 'neuronRemoveList');
end

