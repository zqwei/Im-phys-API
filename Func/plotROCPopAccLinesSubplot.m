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

function plotROCPopAccLinesSubplot(nDataSet, params)

    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);     
    numPlots            = length(nDataSet);
    nBins               = 10;    
    histXout            = 0:0.02:1;
    histFreq            = zeros(length(histXout),length(timePoints) -1);
    
    for nPeriods        = 1:length(timePoints) -1
        nPeriodData     = dataInPeriods(nDataSet, timePoints, nPeriods);     
        areaInt         = zeros(numPlots, 1);
        for nPlot       = 1:numPlots
            nRocTData        = [nPeriodData(nPlot).unit_yes_trial; nPeriodData(nPlot).unit_no_trial];
            nRocOData        = [ones(size(nPeriodData(nPlot).unit_yes_trial)); zeros(size(nPeriodData(nPlot).unit_no_trial))];
            [tp, fp]         = RocFunc(nRocTData, nRocOData, nBins);
            areaInt(nPlot)   = intXY([tp, 1], [fp, 1]);
            areaInt(nPlot)   = max(areaInt(nPlot), 1-areaInt(nPlot));
        end
        histFreq(:,nPeriods) = hist(areaInt, histXout);
        histFreq(:,nPeriods) = histFreq(:,nPeriods)/sum(histFreq(:,nPeriods));
%         histFreq(:,nPeriods) = cumsum(histFreq(:,nPeriods))/sum(histFreq(:,nPeriods));
    end
    
    nTitles = {'pre-sample','sample','delay','response'};
    
    for nPeriods = 2:4
        subplot(1, 3, nPeriods-1)
        hold on
        plot(histXout, histFreq(:, nPeriods),'-', 'linewid', 1.0);
        ylabel('% Accumulated Units')
        xlabel('Area under ROC');
        xlim([0 1])
        ylim([0 1])
        title(nTitles{nPeriods})
        box off
%         set(gca, 'TickDir', 'out')
    end
end


function areaInt       = intXY(vec_x, vec_y)
    
    areaInt            = sum((vec_y(1:end-1)+vec_y(2:end)).*(vec_x(2:end)-vec_x(1:end-1))/2);
    
end

