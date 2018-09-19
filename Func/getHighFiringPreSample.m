%
% getHighFiringPreSample.m
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

function indexNeuron = getHighFiringPreSample(nDataSet, params)

    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);    
    indexNeuron         = false(length(nDataSet), 2);

    for nCell           = 1:length(nDataSet)
        nCellDat        = nDataSet(nCell);
        indexNeuron(nCell, 1) = ttest2(mean(nCellDat.unit_yes_trial(:, 1:timePoints(2))), ...
                                mean(nCellDat.unit_yes_trial(:, timePoints(2)+1:end)), 'tail', 'right');
        indexNeuron(nCell, 2) = ttest2(mean(nCellDat.unit_no_trial(:, 1:timePoints(2))), ...
                                mean(nCellDat.unit_no_trial(:, timePoints(2)+1:end)), 'tail', 'right');
    end

end