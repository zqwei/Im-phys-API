%
% plotHistActivityPop.m
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

function nonActiveNeuronIndex = findNonActiveNeurons(nDataSet, params)
    
    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);                
    yesTrialDataSet     = zeros(length(nDataSet), length(timePoints)-1);
    noTrialDataSet      = zeros(length(nDataSet), length(timePoints)-1);
    tYesTrialDataSet    = false(length(nDataSet), length(timePoints)-1);
    tNoTrialDataSet     = false(length(nDataSet), length(timePoints)-1);
    
    for nPeriods        = 1: length(timePoints) -1        
        nPeriodData     = dataInPeriods(nDataSet, timePoints, nPeriods);
        yesTrialData    = cell2mat(arrayfun(@(x) mean(x.unit_yes_trial), nPeriodData, 'UniformOutput', false));
        tYesTrialData   = cell2mat(arrayfun(@(x) ttest(x.unit_yes_trial(:)), nPeriodData, 'UniformOutput', false));
        tYesTrialData(isnan(tYesTrialData)) = 0;
        yesTrialDataSet(:, nPeriods)  = yesTrialData;
        tYesTrialDataSet(:, nPeriods) = tYesTrialData;        
        noTrialData     = cell2mat(arrayfun(@(x) mean(x.unit_no_trial), nPeriodData, 'UniformOutput', false));
        tNoTrialData    = cell2mat(arrayfun(@(x) ttest(x.unit_no_trial(:)), nPeriodData, 'UniformOutput', false));
        tNoTrialData(isnan(tNoTrialData)) = 0;
        noTrialDataSet(:, nPeriods)   = noTrialData;
        tNoTrialDataSet(:, nPeriods)  = tNoTrialData;
    end    
    
    
    yesNonActiveNeurons  = sum(tYesTrialDataSet,2) == 0;
    noNonActiveNeurons   = sum(tYesTrialDataSet,2) == 0;
    
    nonActiveNeuronIndex = yesNonActiveNeurons & noNonActiveNeurons;
    
end