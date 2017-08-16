%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

maxNumOutlier = 10;


for nData             = 2:6
    load([TempDatDir DataSetList(nData).name '.mat'])
    for nUnit = 1:length(nDataSet)
        nCellDFF = nDataSet(nUnit).unit_yes_trial;
        X = nCellDFF;
        for nTime = 1:size(nCellDFF, 2)
            X(:, nTime) = thomspon_tau(nCellDFF(:, nTime), maxNumOutlier);
        end
        nDataSet(nUnit).unit_yes_trial_removeoutlier = X;
        
        nCellDFF = nDataSet(nUnit).unit_no_trial;
        X = nCellDFF;
        for nTime = 1:size(nCellDFF, 2)
            X(:, nTime) = thomspon_tau(nCellDFF(:, nTime), maxNumOutlier);
        end   
        nDataSet(nUnit).unit_no_trial_removeoutlier = X;
    end
    
    save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet')
   
end

%%
% remove those with above thres
% remove the units with different presample activity

thres          = 0.3;

for nData              = 2:6
    load([TempDatDir DataSetList(nData).name '.mat'])
    params             = DataSetList(nData).params;
    timePeriod         = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    neuronRemoveList   = false(length(nDataSet), 1);
    newDataSet         = rmfield(nDataSet, {'unit_no_trial_removeoutlier', 'unit_yes_trial_removeoutlier'});
    for nUnit          = 1:length(nDataSet)
        nCellData      = nDataSet(nUnit);
        idx_yes_remove = mean(isnan(nCellData.unit_yes_trial_removeoutlier), 2)>thres;
        idx_no_remove  = mean(isnan(nCellData.unit_no_trial_removeoutlier), 2)>thres; 
        if sum(idx_yes_remove)>0 || sum(idx_no_remove)>0 
            if sum(~idx_yes_remove)>20 && sum(~idx_no_remove)>20
                nCellData.unit_yes_trial     = nCellData.unit_yes_trial(~idx_yes_remove, :);
%                 nCellData.unit_yes_trial_raw = nCellData.unit_yes_trial_raw(~idx_yes_remove, :);
                nCellData.unit_yes_trial_index = nCellData.unit_yes_trial_index(~idx_yes_remove);
%                 nCellData.unit_yes_trial_removeoutlier = nCellData.unit_yes_trial_removeoutlier(~idx_yes_remove, :);
                nCellData.unit_no_trial     = nCellData.unit_no_trial(~idx_no_remove, :);
%                 nCellData.unit_no_trial_raw = nCellData.unit_no_trial_raw(~idx_no_remove, :);
                nCellData.unit_no_trial_index = nCellData.unit_no_trial_index(~idx_no_remove);
%                 nCellData.unit_no_trial_removeoutlier = nCellData.unit_no_trial_removeoutlier(~idx_no_remove, :);
                nCellData = rmfield(nCellData, {'unit_no_trial_removeoutlier', 'unit_yes_trial_removeoutlier'});
                if ~ttest2(mean(nCellData.unit_yes_trial(:, timePeriod(1):timePeriod(2))), mean(nCellData.unit_no_trial(:, timePeriod(1):timePeriod(2))))
                    baseline        = nCellData.baseline;
                    rawYes          = (nCellData.unit_yes_trial(:, timePeriod(1):timePeriod(2)) + 1) * baseline;
                    rawNo           = (nCellData.unit_no_trial(:, timePeriod(1):timePeriod(2)) + 1) * baseline;
                    meanPresample   = mean(mean([rawYes; rawNo]));
                    nCellData.unit_yes_trial = (nCellData.unit_yes_trial +1)*baseline/meanPresample-1;
                    nCellData.unit_no_trial  = (nCellData.unit_no_trial +1)*baseline/meanPresample-1;
                    newDataSet(nUnit) = nCellData;
                else
                    neuronRemoveList(nUnit) = true;
                end
            else
                neuronRemoveList(nUnit) = true;
            end
        end  
    end
    sum(~neuronRemoveList)
    newDataSet = newDataSet(~neuronRemoveList);
    nDataSet   = newDataSet;
    save([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'nDataSet', 'neuronRemoveList');
end

