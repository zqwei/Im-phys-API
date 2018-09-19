%
% spikeToImaging.m
%
% This performs the transformation from the spiking data to the Ca++
% imaging data
% 
% Input:
% SpikingData  --- spiking time of a spike train of one cell in one trial  
%
% ----------------------------
% Output:
% CaImaging    --- Tx1 data
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function CaImaging  = spikeTimeToImagingSigmoidNL(spikeTimes, timeSeriesData, params, rMean)
    
    Fm          = params(1);
    Ca0         = params(2);
    n           = params(3);
    tau_decay   = params(4);
    tau_rise    = params(5);
    intNoise    = params(6);
    Ca                        = zeros(length(spikeTimes), length(timeSeriesData));  
    nPreSamplePoints          = 100;
    preSamplePoints           = exprnd(1/rMean, length(spikeTimes), nPreSamplePoints);
    preSamplePoints           = cumsum(preSamplePoints, 2);
    
    for nTrial                = 1:length(spikeTimes)
        % adding n_Pre_Sample_Points
        tSpikeTimes           = spikeTimes{nTrial};
        tPreSamplePoints      = tSpikeTimes(1) - preSamplePoints(nTrial, :);
        tSpikeTimes           = [tPreSamplePoints(end:-1:1), tSpikeTimes'];
        for nSpike            = 1:length(tSpikeTimes)
            Delta_t           = max(timeSeriesData - tSpikeTimes(nSpike),0);
            Ca(nTrial,:)      = Ca(nTrial,:) + exp(-Delta_t/tau_decay).*(1-exp(-Delta_t/tau_rise));
        end
    end
    
    CaIntNoise                = Ca + randn(size(Ca))*intNoise;    
    CaIntNoise(CaIntNoise<0)  = 0;    
    CaImaging                 = Fm./ (1 + exp(-(CaIntNoise-Ca0)*n));