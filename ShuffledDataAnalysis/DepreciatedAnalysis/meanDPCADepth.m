%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;


numComps             = 15;
if ~exist([PlotDir 'CollectedUnitsdPCA'],'dir')
    mkdir([PlotDir 'CollectedUnitsdPCA'])
end

for nData              = 1:length(DataSetList)-1
    load([TempDatDir DataSetList(nData).name '.mat']);
    
    depth                        = [DataSetList(nData).cellinfo(:).depth];    
    depthStart                   = 100;
    depthBin                     = 50;
    depthEnd                     = 900;    
    depth                        = floor((depth-depthStart)/depthBin)*depthBin+depthStart;
    depth(depth>depthEnd)        = depthEnd;
    depth(depth<depthStart)      = depthStart;
    
    uniqueDepth                  = depthStart:depthBin:depthEnd;
    
    
    for nDepth         = 1:length(uniqueDepth)
        if sum(depth == uniqueDepth(nDepth))>50  
            firingRates        = generateDPCAData(nDataSet(depth == uniqueDepth(nDepth)), numTrials);
            firingRatesAverage = nanmean(firingRates, ndims(firingRates));
            trialNum           = ones(size(firingRatesAverage, 1), size(firingRatesAverage, 2))*numTrials;
            optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, ...
                        trialNum, ...
                        'numComps', numComps, ....
                        'combinedParams', combinedParams, ...
                        'numRep', 10, ...  % increase this number to ~10 for better accuracy
                        'display','no');
            [W,V,whichMarg] = dpca(firingRatesAverage, numComps, ...
                        'combinedParams', combinedParams, ...
                        'lambda', optimalLambda);
            explVar = dpcaEV(firingRatesAverage, W, 'combinedParams', combinedParams);

            figure;
            subplot(1, 3, 1)
            plot(1:numComps, [explVar.PCAcomponentVar(1:numComps)', explVar.dPCAcomponentVar'], '-o');
            legend({'PCA', 'dPCA'})
            box off
            xlim([0 numComps+1])
            set(gca,'xTick', 0:5:numComps)
            xlabel('Component index')
            ylabel('% EV')

            subplot(1, 3, 2)
            bar(1:numComps, explVar.PCAmargVar(:, 1:numComps)', 'stack', 'edgecolor', 'none');
            box off
            legend(margNames)
            xlim([0 numComps+1])
            xlabel('Component index')
            set(gca,'xTick', 0:5:numComps)
            ylabel('% EV per PC')

            subplot(1, 3, 3)
            bar(1:numComps, explVar.dPCAmargVar', 'stack', 'edgecolor', 'none');
            box off
            legend(margNames)
            xlim([0 numComps+1])
            xlabel('Component index')
            set(gca,'xTick', 0:5:numComps)
            ylabel('% EV per dPC')

            setPrint(8*3, 6, [PlotDir 'CollectedUnitsdPCA/CollectedUnitsdPCA_' DataSetList(nData).name '_Depth_' num2str(uniqueDepth(nDepth))], 'pdf')

            figure
            bar(1:numComps, explVar.PCAmargVar(:, 1:numComps)', 'stack', 'edgecolor', 'none');
            box off
            legend(margNames)
            xlim([0 numComps+1])
            xlabel('Component index')
            set(gca,'xTick', 0:5:numComps)
            ylabel('% EV per PC')
            setPrint(8, 6, [PlotDir 'CollectedUnitsdPCA/CollectedUnitsPCA_' DataSetList(nData).name '_Depth_' num2str(uniqueDepth(nDepth))], 'pdf')
        end
    end
    
end

close all