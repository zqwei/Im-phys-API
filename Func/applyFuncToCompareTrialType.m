%
% applyFuncToCompareTrialType.m
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

function Results = applyFuncToCompareTrialType (nDataSet, func, funcParams)
    if nargin == 3
        Results = cell2mat(arrayfun(@(x) func(x.unit_yes_trial, x.unit_no_trial, funcParams), nDataSet,'Uniformoutput',false));
    else
        Results = cell2mat(arrayfun(@(x) func(x.unit_yes_trial, x.unit_no_trial), nDataSet,'Uniformoutput',false));
    end