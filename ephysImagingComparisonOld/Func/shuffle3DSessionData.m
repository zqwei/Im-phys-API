%
% shuffle3DSessionData.m
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

function nDataSet = shuffle3DSessionData(nDataSet, totTargets)

    % nDataSet -> numTrial, numUnit, numT
    
    numUnits          = size(nDataSet, 2);
    numYesTrial       = sum(totTargets);
    numNoTrial        = sum(~totTargets);
    
       
    
    for nUnit         = 1:numUnits    
        nUnitYesData  = squeeze(nDataSet(totTargets, nUnit, :));
        nUnitNoData   = squeeze(nDataSet(~totTargets, nUnit, :));
        nDataSet(totTargets, nUnit, :) = nUnitYesData(randperm(numYesTrial) , :);
        nDataSet(~totTargets, nUnit, :) = nUnitNoData(randperm(numNoTrial) , :);
    end
    
end
