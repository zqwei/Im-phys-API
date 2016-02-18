%
% getLogPValueSpikeEpoch.m
%
%
% ----------------------------
% Output:
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%

function logPValue     = getLogPValueTscoreSpikeEpoch(nDataSet, params)

    timePoints         = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    numUnit            = length(nDataSet);
    logPValue          = zeros(numUnit, length(timePoints)-2); % LR
    
    for nUnit          = 1:numUnit
        unit_yes_trial = nDataSet(nUnit).unit_yes_trial;
        unit_no_trial  = nDataSet(nUnit).unit_no_trial;
                                                
        for nPeriod   = 1:length(timePoints)-2
            smoothedYesData = mean(unit_yes_trial(:, timePoints(nPeriod+1):timePoints(nPeriod+2)), 2);
            smoothedNoData  = mean(unit_no_trial(:, timePoints(nPeriod+1):timePoints(nPeriod+2)), 2);
            signData        = sign(mean(smoothedYesData) - mean(smoothedNoData));
            logPValueUnit   = ttest2(smoothedYesData, smoothedNoData) * signData;     
            logPValue(nUnit, nPeriod) = logPValueUnit;           
        end
    end