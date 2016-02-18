%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numFold        = 10;
numComps       = 3;
if ~exist([PlotDir 'CollectedUnitsdPCA'],'dir')
    mkdir([PlotDir 'CollectedUnitsdPCA'])
end

for nData              = [1 3 4]
    
    load([TempDatDir DataSetList(nData).name '.mat']);
    evMat              = zeros(numFold, length(combinedParams), numComps);
    
    figure;
    
    %%% whole population
    firingRates        = generateDPCAData(nDataSet, numTrials);
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));

    pcaX               = firingRatesAverage(:,:);
    firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
    pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));

    % marginalizing
    Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
    % total variance
    totalVar           = sum(sum(pcaX.^2));

    % PCA explained variance
    [~, ~, Wpca] = svd(pcaX');

    PCAmargVar         = zeros(length(combinedParams), length(nDataSet));

    for i=1:length(Xmargs)
        PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar * 100;
    end
    PCAmargVar          = bsxfun(@rdivide, PCAmargVar, sum(PCAmargVar));

    subplot(1, 4, 1);
    plot(1:numComps, PCAmargVar(:, 1:numComps)', '-o', 'linewid', 2)
    box off
    xlim([0.5 3.5])
    ylim([0 1])
    legend(margNames)
%     legend('location','best')
    legend boxoff
    xlabel('Component index')
    ylabel('frac. EV per PC')
    title('Whole population')
    
    
    
    numSubUnit         = 100;
    for nFold          = 1:numFold        
        randIndex          = randperm(length(nDataSet), numSubUnit);    
        firingRates        = generateDPCAData(nDataSet(randIndex), numTrials);
        firingRatesAverage = nanmean(firingRates, ndims(firingRates));

        pcaX               = firingRatesAverage(:,:);
        firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
        pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));

        % marginalizing
        Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
        % total variance
        totalVar           = sum(sum(pcaX.^2));

        % PCA explained variance
        [~, ~, Wpca] = svd(pcaX');

        PCAmargVar         = zeros(length(combinedParams), numSubUnit);

        for i=1:length(Xmargs)
            PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar * 100;
        end
        PCAmargVar          = bsxfun(@rdivide, PCAmargVar, sum(PCAmargVar));

        evMat(nFold, :, :)  = PCAmargVar(:, 1:numComps);
    end

    subplot(1, 4, 2);
    errorbar((1:numComps)'*ones(length(combinedParams), 1)', squeeze(mean(evMat))', squeeze(std(evMat))'/sqrt(numFold), '-o', 'linewid', 2)
    box off
    xlim([0.5 3.5])
    ylim([0 1])
%     legend(margNames)
    xlabel('Component index')
    ylabel('frac. EV per PC')
    title(['Subpopulation n=' num2str(numSubUnit)])
    
    
    numSubUnit         = 500;
    for nFold          = 1:numFold        
        randIndex          = randperm(length(nDataSet), numSubUnit);    
        firingRates        = generateDPCAData(nDataSet(randIndex), numTrials);
        firingRatesAverage = nanmean(firingRates, ndims(firingRates));

        pcaX               = firingRatesAverage(:,:);
        firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
        pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));

        % marginalizing
        Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
        % total variance
        totalVar           = sum(sum(pcaX.^2));

        % PCA explained variance
        [~, ~, Wpca] = svd(pcaX');

        PCAmargVar         = zeros(length(combinedParams), numSubUnit);

        for i=1:length(Xmargs)
            PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar * 100;
        end
        PCAmargVar          = bsxfun(@rdivide, PCAmargVar, sum(PCAmargVar));

        evMat(nFold, :, :)  = PCAmargVar(:, 1:numComps);
    end

    subplot(1, 4, 3);
    errorbar((1:numComps)'*ones(length(combinedParams), 1)', squeeze(mean(evMat))', squeeze(std(evMat))'/sqrt(numFold), '-o', 'linewid', 2)
    box off
    xlim([0.5 3.5])
    ylim([0 1])
%     legend(margNames)
    xlabel('Component index')
    ylabel('frac. EV per PC')
    title(['Subpopulation n=' num2str(numSubUnit)])
    
    
    numSubUnit         = 1000;
    for nFold          = 1:numFold        
        randIndex          = randperm(length(nDataSet), numSubUnit);    
        firingRates        = generateDPCAData(nDataSet(randIndex), numTrials);
        firingRatesAverage = nanmean(firingRates, ndims(firingRates));

        pcaX               = firingRatesAverage(:,:);
        firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
        pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));

        % marginalizing
        Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
        % total variance
        totalVar           = sum(sum(pcaX.^2));

        % PCA explained variance
        [~, ~, Wpca] = svd(pcaX');

        PCAmargVar         = zeros(length(combinedParams), numSubUnit);

        for i=1:length(Xmargs)
            PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar * 100;
        end
        PCAmargVar          = bsxfun(@rdivide, PCAmargVar, sum(PCAmargVar));

        evMat(nFold, :, :)  = PCAmargVar(:, 1:numComps);
    end

    subplot(1, 4, 4);
    errorbar((1:numComps)'*ones(length(combinedParams), 1)', squeeze(mean(evMat))', squeeze(std(evMat))'/sqrt(numFold), '-o', 'linewid', 2)
    box off
    xlim([0.5 3.5])
    ylim([0 1])
%     legend(margNames)
    xlabel('Component index')
    ylabel('frac. EV per PC')
    title(['Subpopulation n=' num2str(numSubUnit)])

    setPrint(8*4, 6, [PlotDir 'CollectedUnitsdPCA/CollectedSubUnitsPCASubAve_' DataSetList(nData).name], 'pdf')
    
end

close all