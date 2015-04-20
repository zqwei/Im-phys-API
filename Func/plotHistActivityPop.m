%
% plotHistActivityPop.m
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

function plotHistActivityPop(nDataSet, params, dist, xlabels, barSeries, yAxes_set, xAxes_set)
    
    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);        
    
    m                   = ceil(sqrt(length(timePoints)-1));
    figure;
    
    % barSeries           = 10;
    
    nTitle              = {'Pre-sample', 'Sample', 'Delay', 'Response'};
    
    for nPeriods        = 1: length(timePoints) -1
        
        subplot(m, m, nPeriods)
        hold on
        nPeriodData     = dataInPeriods(nDataSet, timePoints, nPeriods);
        yesTrialData    = cell2mat(arrayfun(@(x) mean(x.unit_yes_trial), nPeriodData, 'UniformOutput', false));
        barHistWithDist(yesTrialData, dist, xlabels, barSeries, 'b', 1); 
        noTrialData     = cell2mat(arrayfun(@(x) mean(x.unit_no_trial), nPeriodData, 'UniformOutput', false));
        barHistWithDist(noTrialData, dist, xlabels, barSeries, 'r', -1); 
        hold off
        title(nTitle{nPeriods});
        xlim([xAxes_set(1) xAxes_set(end)])
        ylim([yAxes_set(1) yAxes_set(2)])
        s = get(gca,'YTick');
        yTickLabels = arrayfun(@num2str, abs(s), 'UniformOutput', false);
        set(gca,'YTickLabel',yTickLabels);
    
    end    
end