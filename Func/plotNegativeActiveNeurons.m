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

function plotNegativeActiveNeurons(nDataSet, nDataSetRaw, params, xlabels)
    
    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);        
    
    m                   = ceil(sqrt(length(timePoints)-1));
    figure;
    
    % barSeries           = 10;
    
    nTitle              = {'Pre-sample', 'Sample', 'Delay', 'Response'};
    
    yesTrialDataSet     = zeros(length(nDataSet), length(timePoints)-1);
    noTrialDataSet      = zeros(length(nDataSet), length(timePoints)-1);
    tYesTrialDataSet    = false(length(nDataSet), length(timePoints)-1);
    tNoTrialDataSet     = false(length(nDataSet), length(timePoints)-1);
    
    for nPeriods        = 1: length(timePoints) -1        
        nPeriodData     = dataInPeriods(nDataSet, timePoints, nPeriods);
        yesTrialData    = cell2mat(arrayfun(@(x) mean(x.unit_yes_trial), nPeriodData, 'UniformOutput', false));
        tYesTrialData   = cell2mat(arrayfun(@(x) ttest(x.unit_yes_trial(:),0, 'tail','left'), nPeriodData, 'UniformOutput', false));
        tYesTrialData(isnan(tYesTrialData)) = 0;
        yesTrialDataSet(:, nPeriods)  = yesTrialData;
        tYesTrialDataSet(:, nPeriods) = tYesTrialData;        
        noTrialData     = cell2mat(arrayfun(@(x) mean(x.unit_no_trial), nPeriodData, 'UniformOutput', false));
        tNoTrialData    = cell2mat(arrayfun(@(x) ttest(x.unit_no_trial(:),0, 'tail','left'), nPeriodData, 'UniformOutput', false));
        tNoTrialData(isnan(tNoTrialData)) = 0;
        noTrialDataSet(:, nPeriods)   = noTrialData;
        tNoTrialDataSet(:, nPeriods)  = tNoTrialData;
    end    
    
    
    yesNonActiveNeurons = sum(tYesTrialDataSet,2) == 0;
    noNonActiveNeurons  = sum(tYesTrialDataSet,2) == 0;
    
    indexExampleNeurons = find(~yesNonActiveNeurons & ~noNonActiveNeurons);

%     yesNegActiveNeurons = sum(tYesTrialDataSet,2) >= 3;
%     noNegActiveNeurons  = sum(tYesTrialDataSet,2) >= 3;
%     
%     indexExampleNeurons = find(yesNegActiveNeurons & noNegActiveNeurons);


    numExampleNeurons   = 9;
    
    for nNeuron         = 1:numExampleNeurons        
        subplot(numExampleNeurons, 2, nNeuron*2 - 1);
        plotMeanActivityTrace (nDataSet, indexExampleNeurons(nNeuron), params, xlabels, 'Time (s)') 
        subplot(numExampleNeurons, 2, nNeuron*2);
        plotMeanActivityTrace (nDataSetRaw, indexExampleNeurons(nNeuron), params, xlabels, 'Time (s)')        
    end
   
    
end


function plotMeanActivityTrace (nDataSet, nNeuron, params, ylabels, xlabels)
    mean_yes      = mean(nDataSet(nNeuron).unit_yes_trial);
    mean_no       = mean(nDataSet(nNeuron).unit_no_trial);

    var_yes       = sem(nDataSet(nNeuron).unit_yes_trial);
    var_no        = sem(nDataSet(nNeuron).unit_no_trial);
    

    hold on;
    shadedErrorBar(params.timeSeries, mean_yes, var_yes, {'-b', 'linewid', 1.0}, 0.5);
    shadedErrorBar(params.timeSeries, mean_no, var_no, {'-r', 'linewid', 1.0}, 0.5);
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    % ylim([yAxes_set(1) yAxes_set(2)]);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off; 
    xlabel(xlabels);
    ylabel(ylabels)
    
end