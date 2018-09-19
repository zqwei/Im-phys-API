%
% plotROCExampleUnits.m
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

function plotPVALUEPopLines(nDataSet, params)

    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);     
    numPlots            = length(nDataSet);
    % barColor            = {'b', 'r', 'k', 'g'};
    m                   = ceil(sqrt(length(timePoints)-1));
    nBins               = 10;
    figure;
    
    histXout            = 0:0.2:20;
    histFreq            = zeros(length(histXout),length(timePoints) -1);
    
    for nPeriods        = length(timePoints) -1 : -1 : 1
        nPeriodData     = dataInPeriods(nDataSet, timePoints, nPeriods);     
        areaInt         = zeros(numPlots, 1);
        for nPlot       = 1:numPlots
%             nRocTData        = [nPeriodData(nPlot).unit_yes_trial; nPeriodData(nPlot).unit_no_trial];
%             nRocOData        = [ones(size(nPeriodData(nPlot).unit_yes_trial)); zeros(size(nPeriodData(nPlot).unit_no_trial))];
%             [tp, fp]         = RocFunc(nRocTData, nRocOData, nBins);
%             areaInt(nPlot)   = intXY([tp, 1], [fp, 1]);
%             % areaInt(nPlot)   = max(areaInt(nPlot), 1-areaInt(nPlot));
        [~, areaInt(nPlot)]  = ttest2(nPeriodData(nPlot).unit_yes_trial, nPeriodData(nPlot).unit_no_trial);
        end
        % barHistWithDist(areaInt, 'Normal', 'ROC', 0:0.05:1, 'k');
        histFreq(:,nPeriods) = hist(-log(areaInt), histXout);
        histFreq(:,nPeriods) = histFreq(:,nPeriods)/sum(histFreq(:,nPeriods));
    end
    
    plot(histXout, histFreq,'-');
    ylabel('% Units')
    xlabel('-log(p)');
    xlim([0 20])
    ylim([0 0.2])
    legend({'pre-sample','sample','delay','response'},'location','northeast')
    legend('boxoff')
    box off
end


% function areaInt       = intXY(vec_x, vec_y)
%     
%     areaInt            = sum((vec_y(1:end-1)+vec_y(2:end)).*(vec_x(2:end)-vec_x(1:end-1))/2);
%     
% end

