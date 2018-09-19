addpath('../Func');
setDir;

load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
params        = DataSetList(3).params;

dataSetNames  = {DataSetList(4).name, DataSetList(3).name, DataSetList(10).name};

figure;
subplot(2, 3, 1)
hold on
dat = [];
for nData = 1:3
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 1000);
    num_  = sum(tmp(:,1));
    dat   = [dat; [tmp(2,:)'/num_', ones(1000,1)*nData]];
end
violinplot(dat(:, 1), dat(:, 2), 'ShowData', false, 'ViolinAlpha', 1, 'EdgeColor', [1 1 1]);
plot(1:3, [performanceMat.mono], 'k+', 'linewid', 1)
plot(1:3, [refEphys.mono, refEphys.mono, refEphysFast.mono], 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 2)
hold on
dat = [];
for nData = 1:3
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 1000);
    num_  = sum(tmp(:,1));
    dat   = [dat; [tmp(3,:)'/num_', ones(1000,1)*nData]];
end
violinplot(dat(:, 1), dat(:, 2), 'ShowData', false, 'ViolinAlpha', 1, 'EdgeColor', [1 1 1]);
plot(1:3, [performanceMat.multi], 'k+', 'linewid', 1)
plot(1:3, [refEphys.multi, refEphys.multi, refEphysFast.multi], 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 3)
hold on
dat = [];
for nData = 1:3
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp   = [analysisMat.peakiness];
    dat   = [dat; [tmp', ones(1000,1)*nData]];
end
violinplot(dat(:, 1), dat(:, 2), 'ShowData', false, 'ViolinAlpha', 1, 'EdgeColor', [1 1 1]);
plot(1:3, [performanceMat.peak], 'k+', 'linewid', 1)
plot(1:3, [refEphys.peak, refEphys.peak, refEphysFast.peak], 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 4)
hold on
dat = [];
for nData = 1:3
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp         = nan(1000, 1);
    for n_      = 1:1000
        PCAVar  = analysisMat(n_).PCAVar;
        tmp(n_) = 1 - PCAVar(1, 2)/sum(PCAVar(1, :));
    end
    dat   = [dat; [tmp, ones(1000,1)*nData]];
    refPCAImage(nData) = 1 - performanceMat(nData).pca(2)/sum(performanceMat(nData).pca);
end
refPCAEphys  = 1 - refEphys.pca(2)/sum(refEphys.pca);
refPCAEphys_ = 1 - refEphysFast.pca(2)/sum(refEphysFast.pca);
violinplot(dat(:, 1), dat(:, 2), 'ShowData', false, 'ViolinAlpha', 1, 'EdgeColor', [1 1 1]);
plot(1:3, refPCAImage, 'k+', 'linewid', 1)
plot(1:3, [refPCAEphys, refPCAEphys, refPCAEphys_], 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 5)
hold on
dat = [];
for nData = 1:3
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp         = nan(1000, 1);
    for n_      = 1:1000
        tmp(n_) = mean(mean(analysisMat(n_).decodability(:, 8:26)));
    end
    dat   = [dat; [tmp, ones(1000,1)*nData]];
    refPCAImage(nData) = 1 - performanceMat(nData).pca(2)/sum(performanceMat(nData).pca);
end
violinplot(dat(:, 1), dat(:, 2), 'ShowData', false, 'ViolinAlpha', 1, 'EdgeColor', [1 1 1]);
plot(1:3, mean([performanceMat.ldaS]), 'k+', 'linewid', 1)
plot(1:3, mean([refEphys.ldaS, refEphys.ldaS, refEphysFast.ldaS]), 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 6)
hold on
dat = [];
for nData = 1:3
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp         = nan(1000, 1);
    for n_      = 1:1000
        tmp(n_) = mean(mean(analysisMat(n_).decodability(:, 26:47)));
    end
    dat   = [dat; [tmp, ones(1000,1)*nData]];
    refPCAImage(nData) = 1 - performanceMat(nData).pca(2)/sum(performanceMat(nData).pca);
end
violinplot(dat(:, 1), dat(:, 2), 'ShowData', false, 'ViolinAlpha', 1, 'EdgeColor', [1 1 1]);
plot(1:3, mean([performanceMat.ldaD]), 'k+', 'linewid', 1)
plot(1:3, mean([refEphys.ldaD, refEphys.ldaD, refEphysFast.ldaD]), 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')

    
setPrint(8*3, 6*2, 'C2S', 'pdf')
