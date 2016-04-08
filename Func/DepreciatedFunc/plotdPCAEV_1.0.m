%
% plotdPCAEV.m
% 
%
% Spiking dataset
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

function plotdPCAEV (nDataSet, numTrials, numComps, combinedParams, matFileName1, matFileName2, figFileName1, figFileName2)

%     firingRates        = generateDPCAData(nDataSet, numTrials);
%     firingRatesAverage = nanmean(firingRates, ndims(firingRates));
%     trialNum           = ones(size(firingRatesAverage, 1), size(firingRatesAverage, 2))*numTrials;
    margNames          = {'Stim', 'Time', 'Inter'};
%     
%     optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, ...
%                 trialNum, ...
%                 'numComps', numComps, ....
%                 'combinedParams', combinedParams, ...
%                 'numRep', 10, ...  % increase this number to ~10 for better accuracy
%                 'filename', matFileName1,...
%                 'display','no');
%     [W,V,whichMarg] = dpca(firingRatesAverage, numComps, ...
%                 'combinedParams', combinedParams, ...
%                 'lambda', optimalLambda);
%     explVar = dpcaEV(firingRatesAverage, W, 'combinedParams', combinedParams);
%     
%     save(matFileName2, 'optimalLambda', 'W', 'V', 'explVar')
    load(matFileName2, 'optimalLambda', 'explVar')
        
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
    
    setPrint(8*3, 6, figFileName1, 'pdf')
    
    figure
    bar(1:numComps, explVar.PCAmargVar(:, 1:numComps)', 'stack', 'edgecolor', 'none');
    box off
    legend(margNames)
    xlim([0 numComps+1])
    xlabel('Component index')
    set(gca,'xTick', 0:5:numComps)
    ylabel('% EV per PC')
    setPrint(8, 6, figFileName2, 'pdf')
    
end

