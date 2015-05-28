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

function plotROCExampleUnits (nDataSet, sampleSeq, params)

    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);        
    numPlots            = length(sampleSeq);    
    m                   = ceil(sqrt(numPlots));    
%     barColorAll         = french(128, 2);
    nBins               = 10;
    figure;
    
    for nPeriods        = 1:length(timePoints) -1
        nPeriodData     = dataInPeriods(nDataSet(sampleSeq), timePoints, nPeriods);
        
        for nPlot       = 1:numPlots
            subplot(m, m, nPlot)
            hold on;
            plot([0 1],[0 1],'--k')
            nRocTData   = [nPeriodData(nPlot).unit_yes_trial; nPeriodData(nPlot).unit_no_trial];
            nRocOData   = [ones(size(nPeriodData(nPlot).unit_yes_trial)); zeros(size(nPeriodData(nPlot).unit_no_trial))];
            [tp, fp]    = RocFunc(nRocTData, nRocOData, nBins);
            plot(tp, fp);
            hold off;
%             if mean(nPeriodData(nPlot).unit_yes_trial) >= mean(nPeriodData(nPlot).unit_no_trial)
%                 barColor       = barColorAll(128 - ((nPeriods-1)*10+1),:);
%                 plot(tp, fp, 'Color', barColor, 'linewid', 1);
%             else
%                 barColor       = barColorAll(((nPeriods-1)*10+1),:);
%                 plot(fp, tp, 'Color', barColor, 'linewid', 1);
%             end          
%             hold off;
        end
        
    end

end

