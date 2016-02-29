%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population PCA variance and EV of trial type over time
%
% fixed number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListSimultaneous.mat']);

numComps = 15;


if ~exist([PlotDir 'SimultaneousUnitsPCA'],'dir')
    mkdir([PlotDir 'SimultaneousUnitsPCA'])
end

for nData              = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat']);
    coeffPCs           = cell(length(nDataSet), 1);
    pcaVars            = cell(length(nDataSet), 1);
    evTrialTypes       = cell(length(nDataSet), 1);
    
    for nSession      = 1:length(nDataSet)
        nYesDataSet   = nDataSet(nSession).unit_yes_trial;
        numYesTrial   = size(nYesDataSet, 2);
        nNoDataSet    = nDataSet(nSession).unit_no_trial;
        numNoTrial    = size(nNoDataSet, 2);
        nSessionData  = [permute(nYesDataSet, [2 1 3]); permute(nNoDataSet, [2 1 3])];
        nSessionTrial = [true(numYesTrial, 1); false(numNoTrial, 1)];
        numT          = size(nSessionData, 3);
        constRatio    = (numYesTrial + numYesTrial^2/numNoTrial)/(numYesTrial+numNoTrial-1);
        
        pcaVar             = nan(numComps, numT);
        evTrialType        = nan(numComps, numT);
        coeffPC            = nan(size(nSessionData, 2), numT);
        
        for nTime          = 1:numT
            nFiringRates       = squeeze(nSessionData(:, :, nTime));
            [~, score,latent,~,explained,~] = pca(nFiringRates, 'NumComponents', min(numComps, size(nSessionData, 2)));
            coeffPC(:, nTime) = pca(nFiringRates,'numComponents',1);
            coeffPC(:, nTime) = coeffPC(:, nTime)/ norm(coeffPC(:, nTime));
            totVar            = sum(latent);
            pcaVar(1:min(numComps, size(nSessionData, 2)), nTime)  = explained(1:min(numComps, size(nSessionData, 2)))/sum(explained);
            evTrialType(1:min(numComps, size(nSessionData, 2)), nTime) = constRatio * mean(score(nSessionTrial, :)).^2/totVar;
        end
        
        coeffPCs{nSession}    = coeffPC;
        pcaVars{nSession}     = pcaVar;
        evTrialTypes{nSession}= evTrialType;
        
        
        figure;
        subplot(1, 2, 1)
        hold on
        imagesc(DataSetList(nData).params.timeSeries, 1:numComps, pcaVar)
        axis xy;
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([1 numComps])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        xlabel('Time (s)');
        ylabel('Component index');
        title('PC Variance')
        colorbar

        subplot(1, 2, 2)
        hold on
        imagesc(DataSetList(nData).params.timeSeries, 1:numComps, evTrialType)
        axis xy;
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([1 numComps])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        xlabel('Time (s)');
        ylabel('Component index');
        title('PC Trial Type EV')
        colorbar
        
        suptitle(['# units: ' num2str(size(nSessionData, 2)) '; # trials: ' num2str(size(nSessionData, 1))])

        setPrint(12*2, 9, [PlotDir 'SimultaneousUnitsPCA/SimultaneousUnitsPCA_' DataSetList(nData).name  '_Session_' num2str(nSession, '%02d')], 'pdf')
        
    end
    
    
    
    save([TempDatDir 'SimultaneousPCA_' DataSetList(nData).name '.mat'], 'coeffPCs', 'pcaVars', 'evTrialTypes')
end

close all
