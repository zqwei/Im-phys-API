%
% applyFuncToTrialType.m
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

function Results = applyFuncToTrialType (nDataSet, trialType, func)
    
    if strcmp( trialType, 'yes')
        Results = cell2mat(arrayfun(@(x) func(x.unit_yes_trial), nDataSet,'Uniformoutput',false));
    elseif strcmp(trialType, 'no')
        Results = cell2mat(arrayfun(@(x) func(x.unit_no_trial), nDataSet,'Uniformoutput',false));
    else
        disp ('Input error for trial type, please indicate yes/no');
        Results = [];
    end