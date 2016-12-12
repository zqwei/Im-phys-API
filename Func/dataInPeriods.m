%
% dataInPeriods.m
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


function nPeriodData = dataInPeriods(nDataSet, timePoints, nPeriods)
    
    nPeriodData      =  repmat(struct('unit_yes_trial', 1, 'unit_no_trial', 1),length(nDataSet), 1);    
    timeSelected     = timePoints(nPeriods):timePoints(nPeriods+1);    
    
    nData            = arrayfun(@(x) mean(x.unit_yes_trial(:,timeSelected), 2), nDataSet,'UniformOutput', false);  
    [nPeriodData(:).unit_yes_trial] = nData{:};
    
    nData            = arrayfun(@(x) mean(x.unit_no_trial(:,timeSelected), 2), nDataSet,'UniformOutput', false);  
    [nPeriodData(:).unit_no_trial] = nData{:};
    