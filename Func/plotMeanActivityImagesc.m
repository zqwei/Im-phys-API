%
% plotMeanActivityImagesc.m
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

function plotMeanActivityImagesc (nDataSet, params)
    
    figure;
    plotMeanActivityImagescTrialType (nDataSet, params, 'yes');
    
    figure;
    plotMeanActivityImagescTrialType (nDataSet, params, 'no');
    
end

function plotMeanActivityImagescTrialType (nDataSet, params, trialType)

    numUnits = length(nDataSet);
    hold on;
    imagesc(params.timeSeries, 1:numUnits, applyFuncToTrialType(nDataSet, trialType, @mean))
    axis xy
    gridxy ([params.polein, params.poleout, 0],[], 'Color','w','Linestyle','--','linewid', 1.0)
    hold off;
    colorbar
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    ylim([1 numUnits]);    
end