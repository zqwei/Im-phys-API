% 
% generateDPCAData.m
% 
% version 1.0
%
% Comparison list
%
% Output:
% SpikeDataSet     --- yDim x nStim x nDecision (usually ignored) x T x N
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function nDataSet     = generateNonLinearDataFromLinear(spikeDataSet, params)

    nDataSet           = spikeDataSet;
    timePoints         = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    constFMean         = 1.0;
    intNoise           = params.intNoise;
    constNoise         = params.extNoise;
    
    for nData          = 1:length(spikeDataSet)
        paramsSet                                = [params.Fm(nData), params.K(nData), params.n(nData), params.tau_d(nData), params.tau_r(nData), intNoise];
        nDataSet(nData).unit_yes_trial           = CaToFluro(nDataSet(nData).unit_yes_trial_linear, paramsSet);
        nDataSet(nData).unit_no_trial            = CaToFluro(nDataSet(nData).unit_no_trial_linear, paramsSet);
        allTrial_correct                         = [nDataSet(nData).unit_yes_trial; nDataSet(nData).unit_no_trial];
        fMean                                    = max(mean(allTrial_correct(:,timePoints(1):timePoints(2))));
        nDataSet(nData).unit_yes_trial           = (nDataSet(nData).unit_yes_trial - fMean);%/(fMean+constFMean);
        nDataSet(nData).unit_no_trial            = (nDataSet(nData).unit_no_trial - fMean);%/(fMean+constFMean);
        nDataSet(nData).unit_yes_trial           = nDataSet(nData).unit_yes_trial + randn(size(nDataSet(nData).unit_yes_trial))*constNoise;
        nDataSet(nData).unit_no_trial            = nDataSet(nData).unit_no_trial  + randn(size(nDataSet(nData).unit_no_trial))*constNoise;
        nDataSet(nData).unit_yes_trial(:,timePoints(1):timePoints(2)) = 0;
        nDataSet(nData).unit_no_trial(:,timePoints(1):timePoints(2))  = 0;
    end
end

function nTrialFluro = CaToFluro(nTrialCa, params)
    Fm          = params(1);
    K           = params(2);
    n           = params(3);
%     tau_decay   = params(4);
%     tau_rise    = params(5);
    intNoise    = params(6);
    nTrialCa                = nTrialCa + randn(size(nTrialCa))*intNoise;    
%     nTrialCa(nTrialCa<0)    = 0; 
%     nTrialCa(nTrialCa<10)   = 10;
%     nTrialCa(nTrialCa>30)   = 30;
%     nTrialFluro             = Fm* nTrialCa.^n ./ (K^n + nTrialCa.^n);
    nTrialFluro  = nTrialCa;
end