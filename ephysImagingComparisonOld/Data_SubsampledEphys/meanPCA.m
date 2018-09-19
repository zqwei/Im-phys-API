%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numFold        = 30;


numComps       = 3;
cmap = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250];

load([TempDatDir DataSetList(1).name '.mat'])
load('caData.mat')


firingRates        = generateDPCAData(nDataSet, numTrials);
firingRatesAverage = nanmean(firingRates, ndims(firingRates));
pcaX               = firingRatesAverage(:,:);
firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
frThress           = [0 1 4 10];
figure
for nThres         = 1:4 % spike count in this case
    frThres        = frThress(nThres);
    load(['validMat_' num2str(frThres, '%02d')], 'validMat')
    pcaBar    = nan(numFold, length(margNames), numComps);
    for nFold = 1:numFold
        nFR   = firingRatesAverage(validMat(:, nFold), :, :);
        pcaNX = pcaX(validMat(:, nFold), :);
        Xmargs             = dpca_marginalize(nFR, 'combinedParams', combinedParams, 'ifFlat', 'yes');
        totalVar           = sum(sum(pcaNX.^2));
        [~, S, Wpca]       = svd(pcaNX');
        PCAmargVar         = zeros(length(combinedParams), sum(validMat(:, nFold)));
        for i=1:length(Xmargs)
            PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar;
        end
        pcaBar(nFold, :, :) = PCAmargVar(:, 1:numComps)';
    end
    
    hold on
    bar((1:numComps)*2 + nThres*0.4 - 0.4, squeeze(mean(pcaBar)), 0.18,'stacked', 'edgecolor', 'none')

end

bar((1:numComps)*2+1.6, caPCA(:, 1:numComps)', 0.18,'stacked', 'edgecolor', 'none')

box off
% xlim([0.5 numComps+0.5])
% ylim([0 0.4])
xlabel('Component index')
ylabel('frac. EV per PC')
% colormap(cmap(1:3, :))
% set(gca, 'xTick', 1:3)
set(gca, 'TickDir', 'out') 
setPrint(8, 6, 'PCA', 'pdf')