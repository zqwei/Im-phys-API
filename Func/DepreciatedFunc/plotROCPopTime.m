%
% plotROCPopTime.m
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

function plotROCPopTime(nDataSet, params, ROCthres, nColor)

    timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries); 
    numUnits            = length(nDataSet);
    numT                = size(nDataSet(1).unit_yes_trial, 2);

    nBins               = 10;
    
    RocData            = zeros(numUnits, numT);
    
    for nT             = 1:numT
        for nUnit      = 1:numUnits
            nRocTData        = [nDataSet(nUnit).unit_yes_trial(:,nT); nDataSet(nUnit).unit_no_trial(:,nT)];
            nRocOData        = [ones(size(nDataSet(nUnit).unit_yes_trial(:,nT))); zeros(size(nDataSet(nUnit).unit_no_trial(:,nT)))];
            [tp, fp]         = RocFunc(nRocTData, nRocOData, nBins);
            RocData(nUnit, nT)   = intXY([tp, 1], [fp, 1]);
            RocData(nUnit, nT)   = max(RocData(nUnit, nT), 1-RocData(nUnit, nT));
        end
    end
    
    hold on;
    mean_RocData        = mean(RocData>ROCthres);
    plot(params.timeSeries, mean_RocData, '-o', 'linewid', 1.0, 'color', nColor);
    plot(params.timeSeries(timePoints), mean_RocData(timePoints), '*', 'linewid', 1.0, 'color', nColor);
    hold off
    ylim([0 1]);
    xlim([params.timeSeries(1) params.timeSeries(end)]);
end


function areaInt       = intXY(vec_x, vec_y)
    
    areaInt            = sum((vec_y(1:end-1)+vec_y(2:end)).*(vec_x(2:end)-vec_x(1:end-1))/2);
    
end

