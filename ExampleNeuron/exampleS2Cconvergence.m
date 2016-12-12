addpath('../Func');
setDir;

load ([TempDatDir 'DataListShuffle.mat']);
nData = 1; % plot raster and psth
load([TempDatDir DataSetList(nData).name '.mat'])
spikeDataSet = nDataSet;   
params               = DataSetList(nData).params;

sigma                         = 0.15 / params.binsize; % 200 ms
filterLength                  = 11;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse); 

numTrials            = 50;
load([TempDatDir 'Modeled_6s_AAV.mat']);
s2cDataSet           = nDataSet;
s2cFiringRates       = generateDPCADataLinearData(s2cDataSet, numTrials);
s2cFiringRatesMean   = squeeze(mean(s2cFiringRates, 4));
s2cFiringRatesMean   = [squeeze(s2cFiringRatesMean(:, 1, :)), squeeze(s2cFiringRatesMean(:, 2, :))];

load([TempDatDir 'Shuffle_Ca_Slow_Short_Delay_Virus_withOLRemoval.mat']);
caDataSet        = nDataSet;
caFiringRates    = generateDPCAData(caDataSet, numTrials);
caFiringRatesIpsi= caFiringRates;
caFiringRatesIpsi(:, 1, :, :) = caFiringRates(:, 2, :, :);
caFiringRatesIpsi(:, 2, :, :) = caFiringRates(:, 1, :, :);
caFiringRatesMean    = squeeze(mean(caFiringRates, 4));
caFiringRatesIpsiMean= squeeze(mean(caFiringRatesIpsi, 4));
caFiringRatesMean    = [squeeze(caFiringRatesMean(:, 1, :)), squeeze(caFiringRatesMean(:, 2, :))];
caFiringRatesIpsiMean= [squeeze(caFiringRatesIpsiMean(:, 1, :)), squeeze(caFiringRatesIpsiMean(:, 2, :))];

rhoMat               = corr(caFiringRatesMean', s2cFiringRatesMean', 'type', 'Spearman');
rhoMatIpsi           = corr(caFiringRatesIpsiMean', s2cFiringRatesMean', 'type', 'Spearman');
[maxRho, maxIndex]   = max(rhoMat, [], 1);


maxCaIndex           = find(hist(maxIndex, 1:max(maxIndex))>5); %529;

neuronIndex          = find(maxIndex == maxCaIndex(22));
for nCell            = 1:length(neuronIndex)
    nCellid          = neuronIndex(nCell);
    subplot(length(neuronIndex), 2, 2*(nCell-1)+1)
    hold on;
    nUnitData        = spikeDataSet(nCellid).unit_yes_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
    nUnitData        = spikeDataSet(nCellid).unit_no_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('Spikes /s');
    xlabel('Time (s)');
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    set(gca, 'TickDir', 'out')
    subplot(length(neuronIndex), 2, 2*(nCell-1)+2)
    hold on;
    nUnitData        = s2cDataSet(nCellid).unit_yes_trial_linear;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
    nUnitData        = s2cDataSet(nCellid).unit_no_trial_linear;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('df/f');
    xlabel('Time (s)');
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    set(gca, 'TickDir', 'out')
end


% neuronIndex          = [508, 106, 633, 542, 237];
% g = @(p,x) p(1) + p(2)./ (1 + exp((p(3)-x)*p(4)));
% param                = [0, 1, 40, 1; 0, 1, 20, 1; 0, 1, 15, 1; 0, 1, 20, 1; 0, 10, 35, 100]; 
% 
neuronIndex          = [508, 508];
g = @(p,x) p(1) + p(2)./ (1 + exp((p(3)-x)*p(4)));
param                = [0, 20, 40, 1; 0, 2, 10, 1]; 

figure;
for nCell            = 1:length(neuronIndex)
    nCellid          = neuronIndex(nCell);
    
    subplot(length(neuronIndex)+1, 3, 3*(nCell-1)+1)
    hold on;
    nUnitData        = spikeDataSet(nCellid).unit_yes_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
    nUnitData        = spikeDataSet(nCellid).unit_no_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('Spikes /s');
    xlabel('Time (s)');
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    set(gca, 'TickDir', 'out')

    subplot(length(neuronIndex)+1, 3, 3*(nCell-1)+2)
    hold on;
    nUnitData        = g(param(nCell, :), s2cDataSet(nCellid).unit_yes_trial_linear);
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
    nUnitData        = g(param(nCell, :), s2cDataSet(nCellid).unit_no_trial_linear);
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('df/f');
    xlabel('Time (s)');
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    set(gca, 'TickDir', 'out')
    
    
    subplot(length(neuronIndex)+1, 3, 3*(nCell-1)+3)
    hold on;
    nUnitData        = g(param(nCell, :), s2cDataSet(nCellid).unit_yes_trial_linear);
    nUnitData        = imagingToSpike(nUnitData);
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
    nUnitData        = g(param(nCell, :), s2cDataSet(nCellid).unit_no_trial_linear);
    nUnitData        = imagingToSpike(nUnitData);
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('df/f');
    xlabel('Time (s)');
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    set(gca, 'TickDir', 'out')    
end

nCellid          = neuronIndex(nCell);

subplot(length(neuronIndex)+1, 3, 3*nCell+1)
hold on;
nUnitData        = spikeDataSet(nCellid).unit_yes_trial;
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
nUnitData        = spikeDataSet(nCellid).unit_no_trial;
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
ylabel('Spikes /s');
xlabel('Time (s)');
xlim([params.timeSeries(1) params.timeSeries(end)]);
set(gca, 'TickDir', 'out')

subplot(length(neuronIndex)+1, 3, 3*nCell+2)
hold on;
nUnitData        = s2cDataSet(nCellid).unit_yes_trial_linear;
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
nUnitData        = s2cDataSet(nCellid).unit_no_trial_linear;
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
ylabel('df/f');
xlabel('Time (s)');
xlim([params.timeSeries(1) params.timeSeries(end)]);
set(gca, 'TickDir', 'out')


subplot(length(neuronIndex)+1, 3, 3*nCell+3)
hold on;
nUnitData        = s2cDataSet(nCellid).unit_yes_trial_linear;
nUnitData        = imagingToSpike(nUnitData);
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
nUnitData        = s2cDataSet(nCellid).unit_no_trial_linear;
nUnitData        = imagingToSpike(nUnitData);
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
ylabel('df/f');
xlabel('Time (s)');
xlim([params.timeSeries(1) params.timeSeries(end)]);
set(gca, 'TickDir', 'out')

setPrint(8*3, 6*3, 'S2C_C2S')

nCellid = 290;
tau_d = 1.4577;
tau_r = 0.0377;
binsize = 1/14.84;
params.Fm = 1;
params.K  = 1;
params.n  = 1;
params.tau_r = tau_r;
params.tau_d = tau_d;
params.intNoise = 0;
params.extNoise = 0;
params.timeSeries = -6:binsize:5;
nCellDataSet  = getFakeCaImagingData(spikeDataSet(nCellid), params);
figure;
hold on;
nUnitData        = nCellDataSet.unit_yes_trial_linear;
nUnitData        = imagingToSpike(nUnitData)/binsize;
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
nUnitData        = nCellDataSet.unit_no_trial_linear;
nUnitData        = imagingToSpike(nUnitData)/binsize;
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
ylabel('Spikes /s');
xlabel('Time (s)');
xlim([-3.1 2]);
ylim([0 50])
set(gca, 'TickDir', 'out')


nCellid = 290;
tau_d = 1.4577;
tau_r = 0.0377;
binsize = 1/14.84;
params.Fm = 1;
params.K  = 1;
params.n  = 1;
params.tau_r = tau_r;
params.tau_d = tau_d;
params.intNoise = 0;
params.extNoise = 0;
params.timeSeries = -6:binsize:5;
pTime = binsize:binsize:length(params.timeSeries)*binsize;
nCellDataSet  = getFakeCaImagingData(spikeDataSet(nCellid), params);
tau_d = 2.6577;
deconv_filter = (1 - exp(-pTime/tau_r)).* exp(-pTime/tau_d);
C = gallery('circul', [deconv_filter, zeros(1, length(deconv_filter))]);
C = C(1:length(deconv_filter), 1:length(deconv_filter));
C = C';
params.timeSeries = params.timeSeries(1):binsize:params.timeSeries(end);
figure;
hold on;
nUnitData        = nCellDataSet.unit_yes_trial_linear;
for nTrial = 1:size(nUnitData, 1)
    nUnitData(nTrial, :)        = inv(C) * nUnitData(nTrial, :)'/binsize+6; %[deconv(nUnitData(nTrial, :), deconv_filter), zeros(1, length(pTime)-1)];
    % nUnitData(nTrial, :)        = C* nUnitData(nTrial, :)'; % conv(nUnitData(nTrial, 1:end-length(pTime)+1), deconv_filter);
end
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
nUnitData        = nCellDataSet.unit_no_trial_linear;
for nTrial = 1:size(nUnitData, 1)
    nUnitData(nTrial, :)        = inv(C) * nUnitData(nTrial, :)'/binsize+6; %[deconv(nUnitData(nTrial, :), deconv_filter), zeros(1, length(pTime)-1)];
    % nUnitData(nTrial, :)        = C* nUnitData(nTrial, :)';
end
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
ylabel('Spikes /s');
xlabel('Time (s)');
xlim([-3.1 2]);
ylim([0 50])
set(gca, 'TickDir', 'out')

%%
nCellid = 290;
tau_d = 1.4577;
tau_r = 0.0377;
binsize = 1/14.84;
params.Fm = 1;
params.K  = 1;
params.n  = 1;
params.tau_r = tau_r;
params.tau_d = tau_d;
params.intNoise = 0;
params.extNoise = 0;
params.timeSeries = -6:binsize:5;
pTime = binsize:binsize:length(params.timeSeries)*binsize;
nCellDataSet  = getFakeCaImagingData(spikeDataSet(nCellid), params);
yesNoise      = randn(size(nCellDataSet.unit_yes_trial_linear)) * 40;
noNoise       = randn(size(nCellDataSet.unit_no_trial_linear)) * 40;
tau_d = 1.4577;
deconv_filter = (1 - exp(-pTime/tau_r)).* exp(-pTime/tau_d);
C = gallery('circul', [deconv_filter, zeros(1, length(deconv_filter))]);
C = C(1:length(deconv_filter), 1:length(deconv_filter));
C = C';
params.timeSeries = params.timeSeries(1):binsize:params.timeSeries(end);
figure;
hold on;
nUnitData        = nCellDataSet.unit_yes_trial_linear + yesNoise;
for nTrial = 1:size(nUnitData, 1)
    nUnitData(nTrial, :)        = inv(C) * nUnitData(nTrial, :)'/binsize+6; %[deconv(nUnitData(nTrial, :), deconv_filter), zeros(1, length(pTime)-1)];
    % nUnitData(nTrial, :)        = C* nUnitData(nTrial, :)'; % conv(nUnitData(nTrial, 1:end-length(pTime)+1), deconv_filter);
end
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
nUnitData        = nCellDataSet.unit_no_trial_linear + noNoise;
for nTrial = 1:size(nUnitData, 1)
    nUnitData(nTrial, :)        = inv(C) * nUnitData(nTrial, :)'/binsize+6; %[deconv(nUnitData(nTrial, :), deconv_filter), zeros(1, length(pTime)-1)];
    % nUnitData(nTrial, :)        = C* nUnitData(nTrial, :)';
end
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
ylabel('Spikes /s');
xlabel('Time (s)');
xlim([-3.1 2]);
ylim([0 90])
set(gca, 'TickDir', 'out')

figure;
hold on;
nUnitData        = nCellDataSet.unit_yes_trial_linear + yesNoise;
nUnitData        = imagingToSpike(nUnitData)/binsize;
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
nUnitData        = nCellDataSet.unit_no_trial_linear + noNoise;
nUnitData        = imagingToSpike(nUnitData)/binsize;
nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
ylabel('Spikes /s');
xlabel('Time (s)');
xlim([-3.1 2]);
ylim([0 50])
set(gca, 'TickDir', 'out')


%%
% params = DataSetList(nData).params;
% params.Fm = 1;
% params.K  = 1;
% params.n  = 1;
% params.tau_r = tau_r;
% params.tau_d = tau_d;
% params.intNoise = 0;
% params.extNoise = 0;
% params.timeSeries = params.timeSeries(1):binsize/10:5;
% nCellDataSet  = getFakeCaImagingData(spikeDataSet(nCellid), params);
% pTime = binsize/10:binsize/10:length(params.timeSeries)*binsize/10;
% deconv_filter = (1 - exp(-pTime/tau_r)).* exp(-pTime/tau_d);
% C = gallery('circul', [deconv_filter, zeros(1, length(deconv_filter))]);
% C = C(1:length(deconv_filter), 1:length(deconv_filter));
% sigma                         = 0.15 / params.binsize *10; % 200 ms
% filterLength                  = 110;
% filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
% filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
% filterInUse                   = filterInUse / sum (filterInUse); 
% 
% figure;
% hold on;
% nUnitData        = nCellDataSet.unit_yes_trial_linear;
% for nTrial = 1:size(nUnitData, 1)
%     nUnitData(nTrial, :)        = inv(C) * nUnitData(nTrial, :)'; %[deconv(nUnitData(nTrial, :), deconv_filter), zeros(1, length(pTime)-1)];
% end
% nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
% shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
%     std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
%     {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
% nUnitData        = nCellDataSet.unit_no_trial_linear;
% for nTrial = 1:size(nUnitData, 1)
%     nUnitData(nTrial, :)        = inv(C) * nUnitData(nTrial, :)'; %[deconv(nUnitData(nTrial, :), deconv_filter), zeros(1, length(pTime)-1)];
% end
% nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
% shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
%     std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
%     {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
% gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
% hold off;
% ylabel('df/f');
% xlabel('Time (s)');
% xlim([params.timeSeries(1) 2]);
% set(gca, 'TickDir', 'out')