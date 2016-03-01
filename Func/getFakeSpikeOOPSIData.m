% 
% obtain the fake ca++ imaging data using spikes from a list of files
% 
% version 1.0
%
% Comparison list
%
% Output:
% SpikeDataSet     --- yDim x 1 cells (yDims number of neurons) 
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 



function [nDataSetFastOOPSI, nDataSetSMCOOPSI] = getFakeSpikeOOPSIData(spikeDataSet, params)
    
    nDataSetFastOOPSI           = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
                                'unit_yes_trial', 1, 'unit_no_trial', 1),length(spikeDataSet), 1);
                            
    nDataSetSMCOOPSI            = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
                                'unit_yes_trial', 1, 'unit_no_trial', 1),length(spikeDataSet), 1);
                            
    
    % set parameters for run_oopsi                        
    V.fast_iter_max    = length(params.timeSeries);
    V.smc_iter_max     = 1000;
    V.dt               = 1/params.frameRate;
    V.preprocess       = 1; % high-pass filter (increase fitting speed)
    V.T                = length(params.timeSeries);
    
    for nData          = 1:length(spikeDataSet)
        
        disp(num2str(nData));
        
        nDataSetFastOOPSI(nData).sessionIndex     = spikeDataSet(nData).sessionIndex;
        nDataSetFastOOPSI(nData).nUnit            = spikeDataSet(nData).nUnit;
        nDataSetFastOOPSI(nData).depth_in_um      = spikeDataSet(nData).depth_in_um;
        nDataSetFastOOPSI(nData).AP_in_um         = spikeDataSet(nData).AP_in_um;
        nDataSetFastOOPSI(nData).ML_in_um         = spikeDataSet(nData).ML_in_um;
        nDataSetFastOOPSI(nData).cell_type        = spikeDataSet(nData).cell_type;
                
        nDataSetSMCOOPSI(nData).sessionIndex      = spikeDataSet(nData).sessionIndex;
        nDataSetSMCOOPSI(nData).nUnit             = spikeDataSet(nData).nUnit;
        nDataSetSMCOOPSI(nData).depth_in_um       = spikeDataSet(nData).depth_in_um;
        nDataSetSMCOOPSI(nData).AP_in_um          = spikeDataSet(nData).AP_in_um;
        nDataSetSMCOOPSI(nData).ML_in_um          = spikeDataSet(nData).ML_in_um;
        nDataSetSMCOOPSI(nData).cell_type         = spikeDataSet(nData).cell_type;
        
        %% unit_yes_trial
        nDataSetFastOOPSI(nData).unit_yes_trial_index = spikeDataSet(nData).unit_yes_trial_index;
        nDataSetSMCOOPSI(nData).unit_yes_trial_index  = spikeDataSet(nData).unit_yes_trial_index;
        [fastData, smcData]                           = imagingToSpike(spikeDataSet(nData).unit_yes_trial, V);
        nDataSetFastOOPSI(nData).unit_yes_trial       = fastData;
        nDataSetSMCOOPSI(nData).unit_yes_trial        = smcData;
        
        %% unit_no_trial
        nDataSetFastOOPSI(nData).unit_no_trial_index  = spikeDataSet(nData).unit_no_trial_index;
        nDataSetSMCOOPSI(nData).unit_no_trial_index   = spikeDataSet(nData).unit_no_trial_index;
        [fastData, smcData]                           = imagingToSpike(spikeDataSet(nData).unit_no_trial, V);
        nDataSetFastOOPSI(nData).unit_no_trial        = fastData;
        nDataSetSMCOOPSI(nData).unit_no_trial         = smcData;
        
        %% unit_yes_error
        nDataSetFastOOPSI(nData).unit_yes_error_index = spikeDataSet(nData).unit_yes_error_index;
        nDataSetSMCOOPSI(nData).unit_yes_error_index  = spikeDataSet(nData).unit_yes_error_index;
        [fastData, smcData]                           = imagingToSpike(spikeDataSet(nData).unit_yes_error, V);
        nDataSetFastOOPSI(nData).unit_yes_error       = fastData;
        nDataSetSMCOOPSI(nData).unit_yes_error        = smcData;
        
        %% unit_no_trial
        nDataSetFastOOPSI(nData).unit_no_error_index  = spikeDataSet(nData).unit_no_error_index;
        nDataSetSMCOOPSI(nData).unit_no_error_index   = spikeDataSet(nData).unit_no_error_index;
        [fastData, smcData]                           = imagingToSpike(spikeDataSet(nData).unit_no_error, V);
        nDataSetFastOOPSI(nData).unit_no_error        = fastData;
        nDataSetSMCOOPSI(nData).unit_no_error         = smcData;        

    end
