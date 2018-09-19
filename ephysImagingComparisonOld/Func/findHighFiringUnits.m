%
% findHighFiringUnits.m
%
%
% ----------------------------
% Output:
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function selectiveIndex = findHighFiringUnits(nDataSet, params, minFiringRate)
    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    yesTrialDataSet     = zeros(length(nDataSet), length(timePoints)-1);
    noTrialDataSet      = zeros(length(nDataSet), length(timePoints)-1);
    tYesTrialDataSet    = false(length(nDataSet), length(timePoints)-1);
    tNoTrialDataSet     = false(length(nDataSet), length(timePoints)-1);
    
    for nPeriods        = 1: length(timePoints) -1        
        nPeriodData     = dataInPeriods(nDataSet, timePoints, nPeriods);
        yesTrialData    = cell2mat(arrayfun(@(x) mean(x.unit_yes_trial), nPeriodData, 'UniformOutput', false));
        yesTrialDataSet(:, nPeriods)  = yesTrialData;
        noTrialData     = cell2mat(arrayfun(@(x) mean(x.unit_no_trial), nPeriodData, 'UniformOutput', false));
        noTrialDataSet(:, nPeriods)   = noTrialData;
        
        if nPeriods > 1
            tYesTrialData   = cell2mat(arrayfun(@(xIndex) ttest(nPeriodData(xIndex).unit_yes_trial(:), ...
                               yesTrialDataSet(xIndex, 1)), 1:length(nPeriodData), 'UniformOutput', false));
            tYesTrialData(isnan(tYesTrialData)) = 0;
            tYesTrialDataSet(:, nPeriods) = tYesTrialData;                
            tNoTrialData    = cell2mat(arrayfun(@(xIndex) ttest(nPeriodData(xIndex).unit_no_trial(:), ...
                                noTrialDataSet(xIndex, 1)), 1:length(nPeriodData), 'UniformOutput', false));
            tNoTrialData(isnan(tNoTrialData)) = 0;
            tNoTrialDataSet(:, nPeriods)  = tNoTrialData;
        end
    end
    
%     selectiveIndex = (sum(yesTrialDataSet(:,2:end)>minFiringRate,2)>0 | sum(noTrialDataSet(:,2:end)>minFiringRate,2)>0) ...
%                      &  sum(tYesTrialDataSet,2)>0 &  sum(tNoTrialDataSet,2)>0;

    selectiveIndex = sum(yesTrialDataSet(:,2:end)>minFiringRate,2)>0 | sum(noTrialDataSet(:,2:end)>minFiringRate,2)>0;
    
end