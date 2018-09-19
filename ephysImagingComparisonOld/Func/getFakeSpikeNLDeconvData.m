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



function DataSetOOPSI = getFakeSpikeNLDeconvData(spikeDataSet, tau_r, tau_d, nlParams, params)
    
    inv_g              = @(p, x) p(3) - 1/p(4) * log(p(2)./(x-p(1)) - 1);
    per_cent           = 0.02;
    
    DataSetOOPSI       = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
                                'unit_yes_trial', 1, 'unit_no_trial', 1),length(spikeDataSet), 1);
    for nData          = 1:length(spikeDataSet)     
        DataSetOOPSI(nData).sessionIndex     = spikeDataSet(nData).sessionIndex;
        DataSetOOPSI(nData).nUnit            = spikeDataSet(nData).nUnit;
        DataSetOOPSI(nData).depth_in_um      = spikeDataSet(nData).depth_in_um;
        DataSetOOPSI(nData).AP_in_um         = spikeDataSet(nData).AP_in_um;
        DataSetOOPSI(nData).ML_in_um         = spikeDataSet(nData).ML_in_um;
        DataSetOOPSI(nData).cell_type        = spikeDataSet(nData).cell_type;
        
        yesUnitData                          = spikeDataSet(nData).unit_yes_trial;
        noUnitData                           = spikeDataSet(nData).unit_no_trial;
        minData       = min([mean(yesUnitData,1), mean(noUnitData,1)]);
        maxData       = max([mean(yesUnitData,1), mean(noUnitData,1)]);
        randCell      = ceil(rand*length(nlParams));
        % paramInv      = squeeze(nlParams(2, randCell, :));
        paramInv      = squeeze(nlParams(randCell, :));
        paramInv(1)   = minData;
        paramInv(2)   = maxData;
        
        
        
        %% unit_yes_trial
        DataSetOOPSI(nData).unit_yes_trial_index = spikeDataSet(nData).unit_yes_trial_index;
        fastData                                 = mean(spikeDataSet(nData).unit_yes_trial);
        fastData(fastData > paramInv(1) + paramInv(2)*(1-per_cent)) = paramInv(1) + paramInv(2)*(1-per_cent); 
        fastData(fastData < paramInv(1) + paramInv(2)*per_cent)     = paramInv(1) + paramInv(2)*per_cent; 
        fastData                                                    = inv_g(paramInv, fastData);
        fastData                                 = imagingToSpikeDeconv(fastData, tau_r(nData), tau_d(nData), params);
        DataSetOOPSI(nData).unit_yes_trial       = fastData;
        %% unit_no_trial
        DataSetOOPSI(nData).unit_no_trial_index  = spikeDataSet(nData).unit_no_trial_index;
        fastData                                 = mean(spikeDataSet(nData).unit_no_trial);
        fastData(fastData > paramInv(1) + paramInv(2)*(1-per_cent)) = paramInv(1) + paramInv(2)*(1-per_cent); 
        fastData(fastData < paramInv(1) + paramInv(2)*per_cent)     = paramInv(1) + paramInv(2)*per_cent; 
        fastData                                                    = inv_g(paramInv, fastData);
        DataSetOOPSI(nData).unit_no_trial        = fastData;
        %% unit_yes_error
        DataSetOOPSI(nData).unit_yes_error_index = spikeDataSet(nData).unit_yes_error_index;
        fastData                                 = imagingToSpikeDeconv(spikeDataSet(nData).unit_yes_error, tau_r(nData), tau_d(nData), params);
        DataSetOOPSI(nData).unit_yes_error       = fastData;
        %% unit_no_trial
        DataSetOOPSI(nData).unit_no_error_index  = spikeDataSet(nData).unit_no_error_index;
        fastData                                 = imagingToSpikeDeconv(spikeDataSet(nData).unit_no_error, tau_r(nData), tau_d(nData), params);
        DataSetOOPSI(nData).unit_no_error        = fastData;
    end
