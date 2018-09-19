%
% getAnovaLogPSpike.m
% 
%
% Spiking dataset
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

function unitGroup = getAnovaLogPSpikeEpochRewardOnly (logPValueEpoch)

    thresLogP     = -log(0.05);
    numUnit       = size(logPValueEpoch, 1);
    unitGroup     = nan(numUnit, 1);

    gAP           = sum(logPValueEpoch(:, 1, :) > thresLogP, 3) > 0;
    gCE           = sum(logPValueEpoch(:, 3, 3) > thresLogP, 3) > 0;


    unitGroup( gAP & ~gCE) = 1; % group for AP
    unitGroup(~gAP &  gCE) = 2; % group for CE
    unitGroup( gAP &  gCE) = 3; % group for AP + LR + CE 
    unitGroup(~gAP & ~gCE & ~isnan(sum(sum(logPValueEpoch, 3), 2))) = 4; % group for none
    
end

