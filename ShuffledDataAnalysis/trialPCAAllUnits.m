%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population PCA variance and EV of trial type over time
%
% All neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

numTrials      = 5000;
numComps       = 5;
trialType      = [true(numTrials, 1); false(numTrials, 1)];

if ~exist([PlotDir 'CollectedUnitsPCA'],'dir')
    mkdir([PlotDir 'CollectedUnitsPCA'])
end

for nData              = 1:length(DataSetList)-1
%     load([TempDatDir DataSetList(nData).name '.mat']);
%     numUnits           = length(nDataSet);
%     numT               = size(nDataSet(1).unit_yes_trial, 2);
%     pcaVar             = nan(numComps, numT);
%     evTrialType        = nan(numComps, numT);
%     firingRates        = generatePCAData(nDataSet, numTrials);
%     
%     
%     for nTime          = 1:numT    
%         nFiringRates       = squeeze(firingRates(nTime, :, :))';    
%         [~,score,latent,~,explained,~] = pca(nFiringRates, 'NumComponents', numComps);
%         totVar            = sum(latent);
%         pcaVar(:, nTime)  = explained(1:numComps)/sum(explained);
%         evTrialType(:, nTime) = mean(score(1:numTrials, :)).^2/totVar;
%     end
%     
%     save([TempDatDir 'PCATimeAllUnit_' DataSetList(nData).name '.mat'], 'pcaVar', 'evTrialType')
%     
    load([TempDatDir 'PCATimeAllUnit_' DataSetList(nData).name '.mat'], 'pcaVar', 'evTrialType')
    
    figure;
    hold on
    imagesc(DataSetList(nData).params.timeSeries, 1:numComps, pcaVar(1:numComps, :))
    axis xy;
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([1 numComps])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('Component index');
%     title('frac. PC Variance')
    colorbar
    setPrint(8, 6, [PlotDir 'CollectedUnitsPCA/CollectedUnitsPCAFracVar_' DataSetList(nData).name], 'pdf')
    
    figure
    hold on
    imagesc(DataSetList(nData).params.timeSeries, 1:numComps, evTrialType(1:numComps, :))
    axis xy;
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([1 numComps])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('Component index');
%     title('PC Trial Type EV')
    colorbar
    
    setPrint(8, 6, [PlotDir 'CollectedUnitsPCA/CollectedUnitsPCATrialEV_' DataSetList(nData).name], 'pdf')
        
end

close all