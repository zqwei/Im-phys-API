%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population PCA variance and EV of trial type over time
%
% fixed number of neuron with normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListSimultaneous.mat']);
numFold        = 30;

if ~exist([PlotDir '/SimultaneousUnitsDecodability'],'dir')
    mkdir([PlotDir '/SimultaneousUnitsDecodability'])
end

numComps        = 2;

ylimMaxSets     = [0.61, 1.0, 0.4, 0.31];

for nData              = [1 3 4]
    
    ylimMax     = ylimMaxSets(nData);
    
    load([TempDatDir DataSetList(nData).name '.mat']);
    mRow = length(nDataSet);
    
    figure;
    for nSession  = 1:length(nDataSet)
        numYesTrial = length(nDataSet(nSession).unit_yes_trial_index);
        numNoTrial  = length(nDataSet(nSession).unit_no_trial_index);
        totTargets  = [true(numYesTrial, 1); false(numNoTrial, 1)];
        numUnits   = length(nDataSet(nSession).nUnit);
        numTrials  = numYesTrial + numNoTrial;
        rYesTrial  = numYesTrial/numTrials;
        rNoTrial   = numNoTrial/numTrials;
        nSessionData  = [nDataSet(nSession).unit_yes_trial; nDataSet(nSession).unit_no_trial];
        nSessionData  = normalizationDim(nSessionData, 2);  
        numT          = size(nSessionData, 3);
        
        %% compute PC-EV of simultaneous recording data
        simPCAVar     = nan(numComps, numT);
        simEVTrialType= nan(numComps, numT);
        
        for nTime          = 1:numT
            [~,score,latent,~,explained,~] = pca(squeeze(nSessionData(:, :, nTime)), 'NumComponents', numComps);
            totVar            = sum(latent);
            simPCAVar(:, nTime)  = explained(1:numComps)/sum(explained);
            simEVTrialType(:, nTime) = (mean(score(totTargets, :)).^2*rYesTrial + mean(score(~totTargets, :)).^2*rNoTrial)/totVar;
        end
        
        %% compute PC-EV of shuffle recording data
        pcaVar             = nan(numFold, numComps, numT);
        evTrialType        = nan(numFold, numComps, numT);
    
        for nFold             = 1:numFold
            nSessionData  = [nDataSet(nSession).unit_yes_trial; nDataSet(nSession).unit_no_trial];
            nSessionData  = shuffle3DSessionData(nSessionData, totTargets);
            nSessionData  = normalizationDim(nSessionData, 2);
            for nTime          = 1:numT
                [~,score,latent,~,explained,~] = pca(squeeze(nSessionData(:, :, nTime)), 'NumComponents', numComps);
                totVar            = sum(latent);
                pcaVar(nFold, :, nTime)  = explained(1:numComps)/sum(explained);
                evTrialType(nFold, :, nTime) = (mean(score(totTargets, :)).^2*rYesTrial + mean(score(~totTargets, :)).^2*rNoTrial)/totVar;
            end
        end
        
        
        %% plots
        subplot(mRow, mCol, (nSession-1)*4 +1 )
        hold on
        plot(DataSetList(nData).params.timeSeries, simPCAVar(1, :), '-k', 'linewid', 1.0);
        plot(DataSetList(nData).params.timeSeries, simPCAVar(2, :), '-m', 'linewid', 1.0);
        axis xy;
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([0 ylimMax])
        set(gca, 'YTick', [0 roundto(ylimMax,1)])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        xlabel('Time (s)');
        ylabel('frac. PC var.');
        set(gca, 'TickDir', 'out')
        title('Simultaneous')
        
        subplot(mRow, mCol, (nSession-1)*4 +2 )
        hold on
        plot(DataSetList(nData).params.timeSeries, simEVTrialType(1, :), '-k', 'linewid', 1.0);
        plot(DataSetList(nData).params.timeSeries, simEVTrialType(2, :), '-m', 'linewid', 1.0);
        axis xy;
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([0 ylimMax])
        set(gca, 'YTick', [0 roundto(ylimMax,1)])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        xlabel('Time (s)');
        ylabel('frac. EV stim.');
        set(gca, 'TickDir', 'out')
        title(['# units: ' num2str(numUnits)])
        
        subplot(mRow, mCol, (nSession-1)*4 +3 )
        hold on
        shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(pcaVar(:, 1, :), 1)),...
            squeeze(std(pcaVar(:, 1, :), [], 1))/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', 'k'}, 0.5);
        
        shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(pcaVar(:, 2, :), 1)),...
            squeeze(std(pcaVar(:, 2, :), [], 1))/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', 'm'}, 0.5);  
        axis xy;
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([0 ylimMax])
        set(gca, 'YTick', [0 roundto(ylimMax,1)])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        xlabel('Time (s)');
        ylabel('frac. PC var.');
        set(gca, 'TickDir', 'out')
        title('Shuffle')
        
        subplot(mRow, mCol, (nSession-1)*4 +4 )
        hold on
        shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(evTrialType(:, 1, :), 1)),...
                squeeze(std(evTrialType(:, 1, :), [], 1))/sqrt(numFold),...
                {'-', 'linewid', 1.0, 'color', 'k'}, 0.5);

        shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(evTrialType(:, 2, :), 1)),...
                squeeze(std(evTrialType(:, 2, :), [], 1))/sqrt(numFold),...
                {'-', 'linewid', 1.0, 'color', 'm'}, 0.5);  
        axis xy;
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([0 ylimMax])
        set(gca, 'YTick', [0 roundto(ylimMax,1)])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        xlabel('Time (s)');
        ylabel('frac. EV stim.');
        set(gca, 'TickDir', 'out')
        title(['# trials: ' num2str(numTrials)])
        
        
    end
           
    setPrint(8*mCol, 6*mRow, [PlotDir 'SimultaneousUnitsDecodability/SimilarityPCA_' DataSetList(nData).name])
end


margNames = {'PC1', 'PC2'};
cmap = {'k', 'm'};
figure;
hold on
for nColor = 1:length(margNames)
    plot(0, nColor, 's', 'color', cmap{nColor}, 'MarkerFaceColor',cmap{nColor},'MarkerSize', 8)
    text(1, nColor, margNames{nColor})
end
xlim([0 10])
hold off
axis off
setPrint(3, 2, [PlotDir 'SimultaneousUnitsDecodability/SimilarityPCA_Label'])

close all
