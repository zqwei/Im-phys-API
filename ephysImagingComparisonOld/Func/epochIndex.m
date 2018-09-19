%
% plotROCExampleUnits.m
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

function EpochIndex   = epochIndex(params)

    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);        
    EpochIndex          = zeros(length(params.timeSeries), 1);    
    
    for nPeriods        = 1: length(timePoints) -1 
        EpochIndex(timePoints(nPeriods):timePoints(nPeriods+1)-1) = nPeriods;
    end
    
    EpochIndex(end) = nPeriods;

end

