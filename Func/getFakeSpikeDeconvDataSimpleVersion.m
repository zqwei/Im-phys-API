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



function DataSetOOPSI = getFakeSpikeDeconvDataSimpleVersion(spikeDataSet, tau_r, tau_d, params)
    
    DataSetOOPSI           = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
                                'unit_yes_trial', 1, 'unit_no_trial', 1),length(spikeDataSet), 1);
    
    for nData          = 1:length(spikeDataSet)     
        DataSetOOPSI(nData).sessionIndex     = spikeDataSet(nData).sessionIndex;
        DataSetOOPSI(nData).nUnit            = spikeDataSet(nData).nUnit;
        
        DataSetOOPSI(nData).unit_yes_trial_index = spikeDataSet(nData).unit_yes_trial_index;
        fastData                                 = imagingToSpikeDeconv(spikeDataSet(nData).unit_yes_trial, tau_r, tau_d, params);
        DataSetOOPSI(nData).unit_yes_trial       = fastData;

        DataSetOOPSI(nData).unit_no_trial_index  = spikeDataSet(nData).unit_no_trial_index;
        fastData                                 = imagingToSpikeDeconv(spikeDataSet(nData).unit_no_trial, tau_r, tau_d, params);
        DataSetOOPSI(nData).unit_no_trial        = fastData;
    end
    
end