%
% imagingToSpike.m
%
% This performs the transformation from the spiking data to the Ca++
% imaging data
% 
% Input:
% nDataSet  --- nTrial x nTime
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

function contData = imagingToSpike(nDataSet)

    
    contData          = nDataSet;    
    numTrial          = size(nDataSet, 1);
    
    
    for nTrial        = 1:numTrial
%         disp(nTrial)
        [cont, ~]     = conttime_oopsi(nDataSet(nTrial, :));
        contData(nTrial, :) = cont.spk;
    end
