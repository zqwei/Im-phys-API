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

function unitGroup = getAnovaLogPSpikeEpoch (logPValueEpoch)

    thresLogP     = -log(0.01);
    numUnit       = size(logPValueEpoch, 1);
    unitGroup     = nan(numUnit, 1);

    gAP           = sum(logPValueEpoch(:, 1, :) > thresLogP, 3) > 0;
    gLR           = sum(logPValueEpoch(:, 2, :) > thresLogP, 3) > 0;
    gCE           = sum(logPValueEpoch(:, 3, :) > thresLogP, 3) > 0;


    unitGroup( gAP & ~gLR & ~gCE) = 1; % group for AP
    unitGroup(~gAP &  gLR & ~gCE) = 2; % group for LR
    unitGroup(~gAP & ~gLR &  gCE) = 3; % group for CE
    unitGroup( gAP &  gLR & ~gCE) = 4; % group for AP + LR
    unitGroup( gAP & ~gLR &  gCE) = 4;%5; % group for AP + CE
    unitGroup(~gAP &  gLR &  gCE) = 4;%6; % group for LR + CE
    unitGroup( gAP &  gLR &  gCE) = 4;%7; % group for AP + LR + CE 
%     unitGroup(~gAP & ~gLR & ~gCE & ~isnan(sum(sum(logPValueEpoch, 3), 2))) = 5;%8; % group for none
    unitGroup(~gAP & ~gLR & ~gCE) = 5;%8; % group for none
    unitGroup(isnan(sum(sum(logPValueEpoch, 3), 2))) = 6;%8; % group for none
    
end

