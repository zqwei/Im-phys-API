%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similarity of PCA and LDA coefficient vectors as function of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;

if ~exist([PlotDir 'SimultaneousUnitsPCALDACorr'],'dir')
    mkdir([PlotDir 'SimultaneousUnitsPCALDACorr'])
end


for nData             = 1:length(DataSetList)
    load([TempDatDir 'SimultaneousLDA_' DataSetList(nData).name '.mat'], 'coeffVals');
    m                 = ceil(sqrt(length(coeffVals)));
    load([TempDatDir 'SimultaneousUnitStats_' DataSetList(nData).name '.mat'], 'numSimultaneousUnits', 'numSimultaneousTrial');
    [ss, sessionIndex] = sortrows([numSimultaneousUnits, numSimultaneousTrial], [-1 -2]);
    
    figure;
    
    for nSession      = 1:length(coeffVals)
        coeffs       = coeffVals{sessionIndex(nSession)};
        subplot(m, m, nSession);
        hold on;
        imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, coeffs'*coeffs);
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    %     caxis([-1 1]);
        colorbar
        axis xy;
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
        box off;
        hold off;
        xlabel('Time (s)')
        ylabel('Time (s)')
        title(['# units: ' num2str(ss(nSession, 1)) '; # trials: ' num2str(ss(nSession, 2))])
    end
    setPrint(8*m, 6*m, [PlotDir 'SimultaneousUnitsPCALDACorr/SimilarityLDALDA_' DataSetList(nData).name], 'tif')
end

close all
