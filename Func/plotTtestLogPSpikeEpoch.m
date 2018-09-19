%
% plotTtestLogPSpikeEpoch.m
% 
%
% Spiking dataset
%
% ----------------------------
% Output:
% 0: nonselective
% 1: homogenous
% 2: heterogenous
% 
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch)

    numUnit       = size(logPValueEpoch, 1);
    unitGroup     = zeros(numUnit, 1);
    
    for nUnit     = 1:length(unitGroup)
        groupIndex = unique(logPValueEpoch(nUnit, :));
        groupIndex(groupIndex==0) = [];
        unitGroup(nUnit) = length(groupIndex);
    end
    
end

