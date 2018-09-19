%
% plotZScoreImagesc.m
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

function plotZScoreImagesc (nDataSet, params)
    
    h = figure;
    numUnits                = length(nDataSet);
    zScores                 = applyFuncToCompareTrialType(nDataSet, @zScore);
    zScores(isnan(zScores)) = 0;    
    similarityValue         = linkage(zScores, 'complete','correlation');
    [~, ~, similaritySort]  = dendrogram(similarityValue,0);
    close(h)
    figure;
    hold on;
    imagesc(params.timeSeries, 1:numUnits, zScores(similaritySort,:))
    axis xy
    colormap(french(128, 2));
    caxis([-8 8]);
    colorbar
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)    
    hold off;
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    ylim([1 numUnits])
    ylabel('Neuronal index')
    xlabel('Time (s)')
end
