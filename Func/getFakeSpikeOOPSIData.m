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



function DataSetOOPSI = getFakeSpikeOOPSIData(spikeDataSet)
    
    DataSetOOPSI           = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
                                'unit_yes_trial', 1, 'unit_no_trial', 1),length(spikeDataSet), 1);
    for nData          = 1:length(spikeDataSet)     
        DataSetOOPSI(nData).sessionIndex     = spikeDataSet(nData).sessionIndex;
        DataSetOOPSI(nData).nUnit            = spikeDataSet(nData).nUnit;
        DataSetOOPSI(nData).depth_in_um      = spikeDataSet(nData).depth_in_um;
        DataSetOOPSI(nData).AP_in_um         = spikeDataSet(nData).AP_in_um;
        DataSetOOPSI(nData).ML_in_um         = spikeDataSet(nData).ML_in_um;
        DataSetOOPSI(nData).cell_type        = spikeDataSet(nData).cell_type;
        
        %% unit_yes_trial
        DataSetOOPSI(nData).unit_yes_trial_index = spikeDataSet(nData).unit_yes_trial_index;
        fastData                                 = imagingToSpike(spikeDataSet(nData).unit_yes_trial);
        DataSetOOPSI(nData).unit_yes_trial       = fastData;
        %% unit_no_trial
        DataSetOOPSI(nData).unit_no_trial_index  = spikeDataSet(nData).unit_no_trial_index;
        fastData                                 = imagingToSpike(spikeDataSet(nData).unit_no_trial);
        DataSetOOPSI(nData).unit_no_trial        = fastData;
        %% unit_yes_error
        DataSetOOPSI(nData).unit_yes_error_index = spikeDataSet(nData).unit_yes_error_index;
        fastData                                 = imagingToSpike(spikeDataSet(nData).unit_yes_error);
        DataSetOOPSI(nData).unit_yes_error       = fastData;
        %% unit_no_trial
        DataSetOOPSI(nData).unit_no_error_index  = spikeDataSet(nData).unit_no_error_index;
        fastData                                 = imagingToSpike(spikeDataSet(nData).unit_no_error);
        DataSetOOPSI(nData).unit_no_error        = fastData;
    end
