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

function plotROCPopCompLines(nDataSet, params)

    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);     
    numPlots            = length(nDataSet);
    nBins               = 10;
    figure;
    areaInt             = zeros(numPlots, 2);    
    for nPeriods        = 1:2
        nPeriodData     = dataInPeriods(nDataSet, timePoints, nPeriods);     
        for nPlot       = 1:numPlots
            nRocTData        = [nPeriodData(nPlot).unit_yes_trial; nPeriodData(nPlot).unit_no_trial];
            nRocOData        = [ones(size(nPeriodData(nPlot).unit_yes_trial)); zeros(size(nPeriodData(nPlot).unit_no_trial))];
            [tp, fp]         = RocFunc(nRocTData, nRocOData, nBins);
            areaInt(nPlot, nPeriods)   = intXY([tp, 1], [fp, 1]);
            % areaInt(nPlot)   = max(areaInt(nPlot), 1-areaInt(nPlot));
        end
    end
    
    plot(areaInt(:,1), areaInt(:,2),'.');
    hold on;
    plot([0 1],[0 1],'--k')
    hold off;
    xlabel('ROC in presample')
    ylabel('ROC in sample');
    xlim([0 1])
    ylim([0 1])
end


function areaInt       = intXY(vec_x, vec_y)
    
    areaInt            = sum((vec_y(1:end-1)+vec_y(2:end)).*(vec_x(2:end)-vec_x(1:end-1))/2);
    
end

