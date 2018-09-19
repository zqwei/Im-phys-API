%
% ROCPop.m
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

function ROCIndex       = ROCPop(nDataSet, params)

    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);     
    numPlots            = length(nDataSet);
    ROCIndex            = zeros(numPlots, length(timePoints)-1);
    nBins               = 10;
    for nPeriods        = length(timePoints) -1 : -1 : 1
        nPeriodData     = dataInPeriods(nDataSet, timePoints, nPeriods);     
        for nPlot       = 1:numPlots
            nRocTData        = [nPeriodData(nPlot).unit_yes_trial; nPeriodData(nPlot).unit_no_trial];
            nRocOData        = [ones(size(nPeriodData(nPlot).unit_yes_trial)); zeros(size(nPeriodData(nPlot).unit_no_trial))];
            [tp, fp]         = RocFunc(nRocTData, nRocOData, nBins);
            areaInt          = intXY([tp, 1], [fp, 1]);
            ROCIndex(nPlot, nPeriods) = max(areaInt, 1 - areaInt);
        end
    end
    
end


function areaInt       = intXY(vec_x, vec_y)    
    areaInt            = sum((vec_y(1:end-1)+vec_y(2:end)).*(vec_x(2:end)-vec_x(1:end-1))/2);    
end