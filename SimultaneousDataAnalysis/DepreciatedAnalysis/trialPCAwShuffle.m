%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population PCA variance and EV of trial type over time
%
% fixed number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;

load ([TempDatDir 'DataListSimultaneous.mat']);

numFold              = 10;
numComps             = 10;


if ~exist([PlotDir 'SimultaneousShuffleUnitsPCA'],'dir')
    mkdir([PlotDir 'SimultaneousShuffleUnitsPCA'])
end

% for nData              = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat']);
%     load([TempDatDir 'SimultaneousUnitStats_' DataSetList(nData).name '.mat'], 'numSimultaneousUnits', 'numSimultaneousTrial');
% 
%     ssToAnalysis       = numSimultaneousUnits>10 & (numSimultaneousTrial> 2* numSimultaneousUnits);
%     nDataSet           = nDataSet(ssToAnalysis);
%     ss                 = [numSimultaneousUnits(ssToAnalysis), numSimultaneousTrial(ssToAnalysis)];
% 
%     coeffPCs           = cell(length(nDataSet), 1+numFold);
%     pcaVars            = cell(length(nDataSet), 1+numFold);
%     evTrialTypes       = cell(length(nDataSet), 1+numFold);
% 
%     for nSession      = 1:length(nDataSet)
%         nYesDataSet   = nDataSet(nSession).unit_yes_trial;
%         numYesTrial   = size(nYesDataSet, 2);
%         nNoDataSet    = nDataSet(nSession).unit_no_trial;
%         numNoTrial    = size(nNoDataSet, 2);
%         nSessionData  = [permute(nYesDataSet, [2 1 3]); permute(nNoDataSet, [2 1 3])];
%         nSessionTrial = [true(numYesTrial, 1); false(numNoTrial, 1)];
%         numT          = size(nSessionData, 3);
%         constRatio    = (numYesTrial + numYesTrial^2/numNoTrial)/(numYesTrial+numNoTrial-1);
% 
%         pcaVar             = nan(numComps, numT);
%         evTrialType        = nan(numComps, numT);
%         coeffPC            = nan(size(nSessionData, 2), numT);
% 
%         for nTime          = 1:numT
%             nFiringRates       = squeeze(nSessionData(:, :, nTime));
%             [~, score,latent,~,explained,~] = pca(nFiringRates, 'NumComponents', numComps);
%             coeffPC(:, nTime) = pca(nFiringRates,'numComponents',1);
%             coeffPC(:, nTime) = coeffPC(:, nTime)/ norm(coeffPC(:, nTime));
%             totVar                         = sum(latent);
%             pcaVar(1:numComps, nTime)      = explained(1:numComps)/sum(explained);
%             evTrialType(1:numComps, nTime) = constRatio * mean(score(nSessionTrial, :)).^2/totVar;
%         end
% 
%         coeffPCs{nSession, 1}    = coeffPC;
%         pcaVars{nSession, 1}     = pcaVar;
%         evTrialTypes{nSession, 1}= evTrialType;
%         tSessionData             = nSessionData;
%         
%         numNeuron                = size(nSessionData, 2);
% 
%         for nFold                = 1:numFold
%             randYesIndex   = randi(numYesTrial, numYesTrial, numNeuron);
%             randNoIndex    = randi(numNoTrial, numNoTrial, numNeuron) + numYesTrial;
%             randIndex      = [randYesIndex; randNoIndex];
%             for nUnit      = 1:numNeuron
%                 tSessionData(:, nUnit, :)      = nSessionData(randIndex(:, nUnit), nUnit, :);
%             end
% 
%             for nTime          = 1:numT
%                 nFiringRates                   = squeeze(tSessionData(:, :, nTime));
%                 [~, score,latent,~,explained,~] = pca(nFiringRates, 'NumComponents', numComps);
%                 coeffPC(:, nTime)              = pca(nFiringRates,'numComponents',1);
%                 coeffPC(:, nTime)              = coeffPC(:, nTime)/ norm(coeffPC(:, nTime));
%                 totVar                         = sum(latent);
%                 pcaVar(1:numComps, nTime)      = explained(1:numComps)/sum(explained);
%                 evTrialType(1:numComps, nTime) = constRatio * mean(score(nSessionTrial, :)).^2/totVar;
%             end
% 
%             coeffPCs{nSession, 1+nFold}        = coeffPC;
%             pcaVars{nSession, 1+nFold}         = pcaVar;
%             evTrialTypes{nSession, 1+nFold}    = evTrialType;
%         end
% 
%     end
% 
%     save([TempDatDir 'SimultaneousShufflePCA_' DataSetList(nData).name '.mat'], 'coeffPCs', 'pcaVars', 'evTrialTypes', 'ss')
% end


% for nData              = 1:length(DataSetList)
%     load([TempDatDir 'SimultaneousShufflePCA_' DataSetList(nData).name '.mat'], 'evTrialTypes', 'ss');
%     
%     for nSession      = 1:size(evTrialTypes, 1)
%         figure;
%         subplot(1, 2, 1)
%         coeffs       = evTrialTypes{nSession, 1}; 
%         hold on;
%         imagesc(DataSetList(nData).params.timeSeries, 1:numComps, coeffs);
%         xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%         ylim([1 numComps]);
%     %     caxis([-1 1]);
%         colorbar
%         axis xy;
%         gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
%         box off;
%         hold off;
%         xlabel('Time (s)')
%         ylabel('Component index')
%         
%         subplot(1, 2, 2)
%         coeffs       = evTrialTypes{nSession, 2}; 
%         hold on;
%         imagesc(DataSetList(nData).params.timeSeries, 1:numComps, coeffs);
%         xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%         ylim([1 numComps]);
%     %     caxis([-1 1]);
%         colorbar
%         axis xy;
%         gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
%         box off;
%         hold off;
%         xlabel('Time (s)')
%         ylabel('Component index')
%         
%         setPrint(8*2, 6, [PlotDir 'SimultaneousShuffleUnitsPCA/SimultaneousUnitsPCA_' DataSetList(nData).name '_Session_' num2str(nSession)], 'pdf')
%     end
% end



markerSet             = {'o', 's', 'v', '^'};

figure;
hold on
for nData             = 1:length(DataSetList)
    load([TempDatDir 'SimultaneousShufflePCA_' DataSetList(nData).name '.mat'], 'evTrialTypes', 'ss');
    trialVar          = zeros(size(evTrialTypes, 1), 1);
    trialEV           = zeros(size(evTrialTypes, 1), 1);
    numT              = size(evTrialTypes{1, 1}, 2);
    corMat            = zeros(1+numFold, numComps, numT);
    for nSession      = 1:size(evTrialTypes, 1)
        for nFold     = 1:numFold+1
            coeffs                     = evTrialTypes{nSession, nFold}; 
            corMat(nFold, :, :)        = coeffs;  
        end
        
        trialEV(nSession)  = mean(mean((corMat(1, :) - mean(corMat(2:numFold+1, :))).^2, 3), 2)/mean(mean(corMat(1, :, :).^2, 3), 2)*100;
        
        trialVar(nSession) = mean(mean(var(corMat(2:numFold+1, :, :)), 3), 2)/mean(mean(corMat(1, :, :).^2, 3), 2)*100;
    end
    scatter(log10(trialVar), log10(trialEV), ss(:,1), ss(:,2), markerSet{nData}, 'filled')
end
hold off
colormap(jet)
xlabel('% trial by trial var.')
ylabel('% estimation error')
set(gca, 'xTick', 0:1:2, 'xTickLabel', {'1', '10', '100'})
set(gca, 'yTick', 0:1:2, 'yTickLabel', {'1', '10', '100'})
setPrint(8, 6, [PlotDir 'SimultaneousShuffleUnitsPCA/SimultaneousUnitsPCA_' DataSetList(nData).name '_AllSession'], 'pdf')



close all
