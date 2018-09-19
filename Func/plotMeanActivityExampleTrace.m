%
% plotMeanActivityExampleTrace.m
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

function plotMeanActivityExampleTrace (nDataSet, sampleSize, params, ylabels)
    mean_yes      = applyFuncToTrialType(nDataSet, 'yes', @mean);
    mean_no       = applyFuncToTrialType(nDataSet, 'no', @mean);

    var_yes       = applyFuncToTrialType(nDataSet, 'yes', @sem);
    var_no        = applyFuncToTrialType(nDataSet, 'no', @sem);
    
    if length(sampleSize) > 1
        sampleSeq     = sampleSize;
        sampleSize    = length(sampleSeq);
    else
        numUnits      = length(nDataSet);
        sampleSeq     = randperm(numUnits);
        sampleSeq     = sampleSeq(1:sampleSize);
    end
       
    m             = ceil (sqrt(sampleSize));
    
    figure;
    
    for    nPlot  = 1:sampleSize        
        subplot(m, m, nPlot)
        hold on;
        nIndex    = sampleSeq(nPlot);
        shadedErrorBar(params.timeSeries, mean_yes(nIndex,:), var_yes(nIndex,:), {'-b', 'linewid', 1.0}, 0.5);
        shadedErrorBar(params.timeSeries, mean_no(nIndex,:), var_no(nIndex,:), {'-r', 'linewid', 1.0}, 0.5);
        gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
        hold off; 
        xlim([params.timeSeries(1) params.timeSeries(end)]);
        xlabel('Time (s)');
        ylabel(ylabels)
    end
    
end

