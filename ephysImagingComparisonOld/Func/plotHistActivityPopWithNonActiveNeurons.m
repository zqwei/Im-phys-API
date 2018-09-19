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

function plotHistActivityPopWithNonActiveNeurons(nDataSet, params, barSeries, nFactor, xlabels, nDataSetName)
    
    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);        
    
    m                   = ceil(sqrt(length(timePoints)-1));
    figure;
    
    nTitle              = {'Pre-sample', 'Sample', 'Delay', 'Response'};
    
    yesTrialDataSet     = zeros(length(nDataSet), length(timePoints)-1);
    noTrialDataSet      = zeros(length(nDataSet), length(timePoints)-1);
    tYesTrialDataSet    = false(length(nDataSet), length(timePoints)-1);
    tNoTrialDataSet     = false(length(nDataSet), length(timePoints)-1);
    
    for nPeriods        = 1: length(timePoints) -1        
        nPeriodData     = dataInPeriods(nDataSet, timePoints, nPeriods);
        yesTrialData    = cell2mat(arrayfun(@(x) mean(x.unit_yes_trial), nPeriodData, 'UniformOutput', false));
        tYesTrialData   = cell2mat(arrayfun(@(x) ttest(x.unit_yes_trial(:)), nPeriodData, 'UniformOutput', false));
        tYesTrialData(isnan(tYesTrialData)) = 0;
        yesTrialDataSet(:, nPeriods)  = yesTrialData;
        tYesTrialDataSet(:, nPeriods) = tYesTrialData;        
        noTrialData     = cell2mat(arrayfun(@(x) mean(x.unit_no_trial), nPeriodData, 'UniformOutput', false));
        tNoTrialData    = cell2mat(arrayfun(@(x) ttest(x.unit_no_trial(:)), nPeriodData, 'UniformOutput', false));
        tNoTrialData(isnan(tNoTrialData)) = 0;
        noTrialDataSet(:, nPeriods)   = noTrialData;
        tNoTrialDataSet(:, nPeriods)  = tNoTrialData;
    end    
    
    
    yesNonActiveNeurons = sum(tYesTrialDataSet,2) == 0;
    noNonActiveNeurons  = sum(tYesTrialDataSet,2) == 0;
    
    for nPeriods        = 1: length(timePoints) -1
        subplot(m, m, nPeriods)                
        hold on        
        nonHistFreq      = histc(yesTrialDataSet(yesNonActiveNeurons, nPeriods), barSeries);
        nonEpochHistFreq = histc(yesTrialDataSet(~yesNonActiveNeurons & ~tYesTrialDataSet(:, nPeriods), nPeriods), barSeries);
        actHistFreq      = histc(yesTrialDataSet(~yesNonActiveNeurons & tYesTrialDataSet(:, nPeriods), nPeriods), barSeries);
        
        barPlot              = bar(barSeries, [actHistFreq, nonEpochHistFreq, nonHistFreq],'stacked', 'edgeColor','b', 'faceColor', 'b','linewid',0.2);%,'linewid',1
        barPlot(2).FaceColor = [0.5 0.5 1];
        barPlot(3).FaceColor = 'w';
%         barPlot              = bar(barSeries, [actHistFreq]', 'edgeColor','b', 'faceColor', 'b');
        
        nonHistFreq      = histc(noTrialDataSet(noNonActiveNeurons, nPeriods), barSeries);
        nonEpochHistFreq = histc(noTrialDataSet(~noNonActiveNeurons & ~tNoTrialDataSet(:, nPeriods), nPeriods), barSeries);
        actHistFreq      = histc(noTrialDataSet(~noNonActiveNeurons & tNoTrialDataSet(:, nPeriods), nPeriods), barSeries);
        
        barPlot              = bar(barSeries', -[actHistFreq, nonEpochHistFreq, nonHistFreq],'stacked', 'edgeColor','r', 'faceColor', 'r','linewid',0.2);%,'linewid',1
        barPlot(2).FaceColor = [1 0.5 0.5];
        barPlot(3).FaceColor = 'w';
%         barPlot              = bar(barSeries, -[actHistFreq]', 'edgeColor','r', 'faceColor', 'r');
        
        hold off;
        title(nTitle{nPeriods});
        % xlim([xAxes_set(1) xAxes_set(end)])
        % ylim([yAxes_set(1) yAxes_set(2)])
        
        xlim([barSeries(1) barSeries(end)])
        ylim([-length(nDataSet)*nFactor length(nDataSet)*nFactor])
        xlabel(xlabels)
        ylabel('Count (#)')
        s = get(gca,'YTick');
        yTickLabels = arrayfun(@num2str, abs(s), 'UniformOutput', false);
        set(gca,'YTick',s,'YTickLabel',yTickLabels);    
    end    
    
end