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

function CaImaging    = spikeToImagingNonCell(SpikingData, params)

    if nargin         == 1
        params        = [];
    end

    if isfield(params, 'Fm')
        Fm            = params.Fm; % Scalar
    else
        Fm            = 20.0000; % Scalar
    end
    
    if isfield(params, 'Kd')
        Kd            = params.Kd; % Constant
    else
        Kd            = 11.5885; % Constant
    end

    if isfield(params, 'n')
        n             = params.n;  % power
    else
        n             = 2.2152;  % power
    end
    
    if isfield(params, 'tau_decay')
        tau_decay     = params.tau_decay;  % sec
    else
        tau_decay     = 1.5471;  % sec
    end
    
    if isfield(params, 'tau_rise')
        tau_rise      = params.tau_rise;  % sec
    else
        tau_rise      = 0.0670;  % sec
    end
    
    if isfield(params, 'binsize')
        binsize       = params.binsize; % sec
    else
        binsize       = 50/1000; % sec
    end
    
    if isfield(params, 'timeSeries')
        timeSeries    = params.timeSeries;
    else
        timeSeries    = 0 : binsize : 5.3;
    end
    
    
        
%     CaImaging          = zeros(length(timeSeries), 1);
    Ca          = zeros(1, length(timeSeries));
    for nSpike  = 1:length(SpikingData)
        Delta_t = max(timeSeries - SpikingData(nSpike),0);
        Ca      = Ca + exp(-Delta_t/tau_decay).*(1-exp(-Delta_t/tau_rise));
    end
    CaImaging   = Fm* Ca.^n ./ (Kd^n + Ca.^n);
