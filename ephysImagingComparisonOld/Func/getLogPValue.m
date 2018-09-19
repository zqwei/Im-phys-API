%
% getLogPValue.m
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

function logPValue     = getLogPValue(nDataSet)

    numUnit            = length(nDataSet);
    numTime            = size(nDataSet(1).unit_yes_trial, 2);
    logPValue          = zeros(numUnit, 3, numTime); % g1 for AP; g2 for LR; g3 for CE
    
    for nUnit          = 1:numUnit
        numYesTrial    = length(nDataSet(nUnit).unit_yes_trial_index);
        numNoTrial     = length(nDataSet(nUnit).unit_no_trial_index);
        numYesError    = length(nDataSet(nUnit).unit_yes_error_index);
        numNoError     = length(nDataSet(nUnit).unit_no_error_index);
        
        nUnitData      = [  nDataSet(nUnit).unit_yes_trial; ...
                            nDataSet(nUnit).unit_no_trial; ...
                            nDataSet(nUnit).unit_yes_error; ...
                            nDataSet(nUnit).unit_no_error];
                        
        gLR            = [  true(numYesTrial, 1); ...
                            false(numNoTrial, 1); ...
                            true(numYesError, 1); ...
                            false(numNoError, 1)];
        
        gAP            = [  true(numYesTrial, 1); ...
                            false(numNoTrial, 1); ...
                            false(numYesError, 1); ...
                            true(numNoError, 1)];
                        
        gCE            = [  true(numYesTrial, 1); ...
                            true(numNoTrial, 1); ...
                            false(numYesError, 1); ...
                            false(numNoError, 1)]; 
                        
        logPValueUnit  = arrayfun(@(x) -log(anovan(nUnitData(:, x), {gAP, gLR, gCE}, ...
                        'display', 'off')), 1:numTime, 'Uniformoutput',false);     
        
        logPValue(nUnit, :, :) = cell2mat(logPValueUnit);            
    end