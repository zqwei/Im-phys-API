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

function [fastData, smcData] = imagingToSpike(nDataSet, V)

    
    fastData          = nDataSet;
    smcData           = nDataSet;
    
    numTrial          = size(nDataSet, 1);
    
    for nTrial        = 1:numTrial
        [fast, smc]   = runOOPSI(nDataSet(nTrial, :), V);
        n_fast        = fast.n/max(fast.n);
        fastData(nTrial, :) = n_fast / V.dt;
        if ~isempty(smc)
            smcData(nTrial, :)  = smc.E.nbar/V.dt;
        else
            smcData(nTrial, :)  = n_fast / V.dt;
        end
    end
