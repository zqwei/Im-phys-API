%
% plotHistActivityExampleUnits.m
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

function plotHistActivityExampleUnits (nDataSet, sampleSeq, params, dist)
    
    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);    
    
    numPlots            = length(sampleSeq);
    
    m                   = ceil(sqrt(numPlots));
    barSeries           = 10;
    
    barColorAll         = french(128, 2);
    
    figure;
    
    for nPeriods        = length(timePoints) -1 : -1 : 1
        nPeriodData     = dataInPeriods(nDataSet(sampleSeq), timePoints, nPeriods);
        
        for nPlot       = 1:numPlots
            subplot(m, m, nPlot)
            barColor    = barColorAll(128 - ((nPeriods-1)*10+1),:);
            hold on;
            barData     = nPeriodData(nPlot).unit_yes_trial;
            barSign     = 1;
            barHistWithDist(barData(:), dist, barSeries, barColor, barSign); 
            barColor    = barColorAll((nPeriods-1)*10+1,:);
            barSign     = -1;
            barHistWithDist(barData(:), dist, barSeries, barColor, barSign); 
            hold off;
        end
        
    end

    
end