% 
% obtain synthetic spikes using ca++ imaging data from a list of files
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



function DataSetOOPSI = getSyntheticSpikeDeconvDataSimpleVersion(caDataSet, tau_r, tau_d, params)
    
    DataSetOOPSI           = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
                                'unit_yes_trial', 1, 'unit_no_trial', 1),length(caDataSet), 1);
    
    for nData          = 1:length(caDataSet)     
        DataSetOOPSI(nData).sessionIndex     = caDataSet(nData).sessionIndex;
        DataSetOOPSI(nData).nUnit            = caDataSet(nData).nUnit;
        
        DataSetOOPSI(nData).unit_yes_trial_index = caDataSet(nData).unit_yes_trial_index;
        fastData                                 = imagingToSpikeDeconv(caDataSet(nData).unit_yes_trial, tau_r, tau_d, params);
        DataSetOOPSI(nData).unit_yes_trial       = fastData;

        DataSetOOPSI(nData).unit_no_trial_index  = caDataSet(nData).unit_no_trial_index;
        fastData                                 = imagingToSpikeDeconv(caDataSet(nData).unit_no_trial, tau_r, tau_d, params);
        DataSetOOPSI(nData).unit_no_trial        = fastData;
    end
    
end