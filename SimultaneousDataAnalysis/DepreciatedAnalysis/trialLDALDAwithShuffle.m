%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Different ROC
% Same number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;

load ([TempDatDir 'DataListSimultaneous.mat']);

if ~exist([PlotDir '/SimultaneousShuffleUnitsPCALDACorr'],'dir')
    mkdir([PlotDir '/SimultaneousShuffleUnitsPCALDACorr'])
end

numFold              = 10;

% for nData             = 1:length(DataSetList)
%     load([TempDatDir 'SimultaneousShuffleLDA_' DataSetList(nData).name '.mat'], 'coeffVals', 'ss');
%     
%     for nSession      = 1:size(coeffVals, 1)
%         figure;
%         subplot(1, 2, 1)
%         coeffs       = coeffVals{nSession, 1}; 
%         hold on;
%         imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, coeffs'*coeffs);
%         xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%         ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     %     caxis([-1 1]);
%         colorbar
%         axis xy;
%         gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
%         box off;
%         hold off;
%         xlabel('Time (s)')
%         ylabel('Time (s)')
%         
%         subplot(1, 2, 2)
%         coeffs       = coeffVals{nSession, 2}; 
%         hold on;
%         imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, coeffs'*coeffs);
%         xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%         ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     %     caxis([-1 1]);
%         colorbar
%         axis xy;
%         gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
%         box off;
%         hold off;
%         xlabel('Time (s)')
%         ylabel('Time (s)')
%         
%         setPrint(8*2, 6, [PlotDir 'SimultaneousShuffleUnitsPCALDACorr/SimilarityLDALDA_' DataSetList(nData).name '_Session_' num2str(nSession)], 'pdf')
%     end
%     
% end



markerSet             = {'o', 's', 'v', '^'};

figure;
hold on
for nData             = 1:length(DataSetList)
    load([TempDatDir 'SimultaneousShuffleLDA_' DataSetList(nData).name '.mat'], 'coeffVals', 'ss');
    trialVar          = zeros(size(coeffVals, 1), 1);
    trialEV           = zeros(size(coeffVals, 1), 1);
    numT              = size(coeffVals{1, 1}, 2);
    corMat            = zeros(1+numFold, numT, numT);
    for nSession      = 1:size(coeffVals, 1)
        for nFold     = 1:numFold+1
            coeffs                     = coeffVals{nSession, nFold}; 
            corMat(nFold, :, :)        = coeffs'*coeffs;  
        end
        
        trialEV(nSession)  = mean(mean((corMat(1, :) - mean(corMat(2:numFold+1, :))).^2, 3), 2)/mean(mean(corMat(1, :, :).^2, 3), 2)*100;
        
        trialVar(nSession) = mean(mean(var(corMat(2:numFold+1, :, :)), 3), 2)/mean(mean(corMat(1, :, :).^2, 3), 2)*100;
    end
    scatter(trialVar, trialEV, ss(:,1), ss(:,2), markerSet{nData}, 'filled')
end
hold off
colormap(jet)
xlabel('% trial by trial var.')
ylabel('% estimation error')
setPrint(8, 6, [PlotDir 'SimultaneousShuffleUnitsPCALDACorr/SimilarityLDALDA_' DataSetList(nData).name '_AllSession'], 'pdf')

close all;