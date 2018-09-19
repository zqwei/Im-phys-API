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



function nDataSet = getFakeCaImagingDataSigmoid(spikeDataSet, params)
    
    nDataSet           = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
                                'unit_yes_trial', 1, 'unit_no_trial', 1),length(spikeDataSet), 1);
                            
    timePoints         = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    timeSeriesData     = params.timeSeries;
    constFMean         = 0.1;
    intNoise           = params.intNoise;
    constNoise         = params.extNoise;
    
    for nData          = 1:length(spikeDataSet)
        paramsSet                                = [params.Fm(nData), params.Ca0(nData), params.n(nData), params.tau_d(nData), params.tau_r(nData), intNoise];
        nDataSet(nData).sessionIndex             = spikeDataSet(nData).sessionIndex;
        nDataSet(nData).nUnit                    = spikeDataSet(nData).nUnit;
        allTrialFirngAct_correct                 = [spikeDataSet(nData).unit_yes_trial; spikeDataSet(nData).unit_no_trial];
        rMean                                    = mean(mean(allTrialFirngAct_correct(:,timePoints(1):timePoints(2))));
        if rMean == 0; rMean = rand(); end % pre-sample is 0.5 sec, if there is one spike, then the minimal rate is 2 Hz.
        [nDataSet(nData).unit_yes_trial, nDataSet(nData).unit_yes_trial_linear]  = spikeTimeToImagingSigmoid(spikeDataSet(nData).unit_yes_trial_spk_time, timeSeriesData, paramsSet, rMean);
        nDataSet(nData).unit_yes_trial_index     = spikeDataSet(nData).unit_yes_trial_index;
        [nDataSet(nData).unit_no_trial, nDataSet(nData).unit_no_trial_linear] = spikeTimeToImagingSigmoid(spikeDataSet(nData).unit_no_trial_spk_time, timeSeriesData, paramsSet, rMean);
        nDataSet(nData).unit_no_trial_index      = spikeDataSet(nData).unit_no_trial_index;
        allTrial_correct                         = [nDataSet(nData).unit_yes_trial; nDataSet(nData).unit_no_trial];
        fMean                                    = mean(mean(allTrial_correct(:,timePoints(1):timePoints(2))));
        nDataSet(nData).unit_yes_trial           = (nDataSet(nData).unit_yes_trial - fMean)/(fMean+constFMean);
        nDataSet(nData).unit_no_trial            = (nDataSet(nData).unit_no_trial - fMean)/(fMean+constFMean);
        nDataSet(nData).unit_yes_trial           = nDataSet(nData).unit_yes_trial + randn(size(nDataSet(nData).unit_yes_trial))*constNoise;
        nDataSet(nData).unit_no_trial            = nDataSet(nData).unit_no_trial  + randn(size(nDataSet(nData).unit_no_trial))*constNoise;
    
        [nDataSet(nData).unit_yes_error, nDataSet(nData).unit_yes_error_linear] = spikeTimeToImagingSigmoid(spikeDataSet(nData).unit_yes_error_spk_time, timeSeriesData, paramsSet, rMean);
        nDataSet(nData).unit_yes_error_index     = spikeDataSet(nData).unit_yes_error_index;
        [nDataSet(nData).unit_no_error , nDataSet(nData).unit_no_error_linear]  = spikeTimeToImagingSigmoid(spikeDataSet(nData).unit_no_error_spk_time, timeSeriesData, paramsSet, rMean);
        nDataSet(nData).unit_no_error_index      = spikeDataSet(nData).unit_no_error_index;  
        nDataSet(nData).unit_yes_error           = (nDataSet(nData).unit_yes_error - fMean)/(fMean+constFMean);
        nDataSet(nData).unit_no_error            = (nDataSet(nData).unit_no_error - fMean)/(fMean+constFMean);
        nDataSet(nData).unit_yes_error           = nDataSet(nData).unit_yes_error+ randn(size(nDataSet(nData).unit_yes_error))*constNoise;
        nDataSet(nData).unit_no_error            = nDataSet(nData).unit_no_error + randn(size(nDataSet(nData).unit_no_error))*constNoise;

        nDataSet(nData).depth_in_um              = spikeDataSet(nData).depth_in_um;
        nDataSet(nData).AP_in_um                 = spikeDataSet(nData).AP_in_um;
        nDataSet(nData).ML_in_um                 = spikeDataSet(nData).ML_in_um;
        if strcmp(spikeDataSet(nData).cell_type, 'putative_interneuron')
            nDataSet(nData).cell_type                = 0;
        elseif strcmp(spikeDataSet(nData).cell_type, 'putative_pyramidal')
            nDataSet(nData).cell_type                = 1;
        elseif isnumeric(spikeDataSet(nData).cell_type)
            nDataSet(nData).cell_type                = spikeDataSet(nData).cell_type;
        else
            nDataSet(nData).cell_type                = nan;
        end
    end
