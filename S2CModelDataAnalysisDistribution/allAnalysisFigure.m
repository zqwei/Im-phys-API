addpath('../Func');
setDir;

load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
params        = DataSetList(1).params;

dataSetNames  = {'Modeled_6s_AAV_nParaSet_', 'Modeled_GP43_nParaSet_'};


nData         = 1;
load([TempDatDir 'ResultsCompiled_' dataSetNames{nData} '.mat'], 'analysisMat')
X             = reshape([analysisMat.pCa0], 30, 30);
Y             = reshape([analysisMat.pTaud], 30, 30);
Z             = reshape([analysisMat.sizeGroup], 3, 900);
figure;
subplot(1, 3, 1)
h = pcolor(X, Y, reshape(Z(2,:), 30, 30)/720-0.66);
set(h, 'EdgeColor', 'none');
colormap(french(128, 2))
caxis([-0.2 0.2])
colorbar
box off
xlim([0 100])
ylim([0 100])
xlabel('Frac. Nonlinearity')
ylabel('Frac. \tau_d')
set(gca, 'TickDir', 'out')
subplot(1, 3, 2)
h = pcolor(X, Y, reshape(Z(3,:), 30, 30)/720-0.05);
colormap(french(128, 2))
caxis([-0.1 0.1])
set(h, 'EdgeColor', 'none');
colorbar
box off
xlim([0 100])
ylim([0 100])
xlabel('Frac. Nonlinearity')
ylabel('Frac. \tau_d')
set(gca, 'TickDir', 'out')
Z             = [analysisMat.peakiness];
subplot(1, 3, 3)
h = pcolor(X, Y, reshape(Z, 30, 30)-0.49);
colormap(french(128, 2))
caxis([-0.2 0.2])
set(h, 'EdgeColor', 'none');
colorbar
box off
xlim([0 100])
ylim([0 100])
xlabel('Frac. Nonlinearity')
ylabel('Frac. \tau_d')
set(gca, 'TickDir', 'out')
setPrint(8*4, 6, ['s2c_example_' dataSetNames{nData}], 'pdf')