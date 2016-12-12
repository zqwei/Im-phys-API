% 
% obtain the fake dff of spikes from a list of files
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



function nDataSet = getDFFSpike(spikeDataSet, params)
    
    nDataSet           = spikeDataSet;
                            
    timePoints         = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    constFMean         = 0;
    
    for nData          = 1:length(spikeDataSet)
        allTrial_correct                         = [nDataSet(nData).unit_yes_trial; nDataSet(nData).unit_no_trial];
        fMean                                    = mean(mean(allTrial_correct(:,timePoints(1):timePoints(2))));
        if fMean == 0; constFMean = rand(); end % pre-sample is 0.5 sec, if there is one spike, then the minimal rate is 0-1 Hz.        
        nDataSet(nData).unit_yes_trial           = (nDataSet(nData).unit_yes_trial - fMean)/(fMean+constFMean);
        nDataSet(nData).unit_no_trial            = (nDataSet(nData).unit_no_trial - fMean)/(fMean+constFMean);
        nDataSet(nData).unit_yes_error           = (nDataSet(nData).unit_yes_error - fMean)/(fMean+constFMean);
        nDataSet(nData).unit_no_error            = (nDataSet(nData).unit_no_error - fMean)/(fMean+constFMean);
    end
