%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Similarity of PCA and LDA coefficient vectors as function of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_10

addpath('../Func');
setDir;
load ([TempDatDir 'DataList.mat']);
addNoise         = [1 0 0 0 0 0];
numTrials           = 50;
% totTargets          = [true(numTrials/2,1); false(numTrials/2,1)];

if ~exist([PlotDir '/Collected_Units_PCA_LDA'],'dir')
    mkdir([PlotDir '/Collected_Units_PCA_LDA'])
end

nData        = 1;
load([TempDatDir DataSetList(nData).name '.mat'])

noise_amp    = 1e-10;


figure;
subplot(1, 3, 1)
title('Simultaneous Data')
hold on
n_index      = 9;
totTargets   = [true(size(nDataSet3D(n_index).unit_yes_trial,2),1);...
                false(size(nDataSet3D(n_index).unit_no_trial,2),1)];
nSessionData = [permute(nDataSet3D(n_index).unit_yes_trial, [2, 1, 3]);...
                permute(nDataSet3D(n_index).unit_no_trial, [2, 1, 3])];
% nSessionData = nSessionData + randn(size(nSessionData))*noise_amp/sqrt(numTrials)* addNoise(nData);
% nSessionData = normalizationDim(nSessionData, 2);
coeffs       = coeffLDA(nSessionData, totTargets);
% coeffs       = coeffLDA([nSessionData; nSessionData], totTargets);
imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, coeffs'*coeffs);
xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
caxis([-1 1]);
axis xy;
gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
box off;
hold off;
xlabel('Time (s)')
ylabel('Time (s)')


% figure;
% 
% m = ceil(sqrt(size(nSessionData,2)));
% for nNeuron = 1:size(nSessionData,2)
%     subplot(m, m, nNeuron)
%     hold on
%     plot(squeeze(mean(nSessionData(totTargets,nNeuron,:),1)),'-b');
%     plot(squeeze(mean(nSessionData(~totTargets,nNeuron,:),1)),'-r');
%     hold off
% end

subplot(1, 3, 2)
title('Simultaneous Data w/ Shuffle')
hold on
simultaneous_index = [nDataSet(:).sessionIndex] == nDataSet3D(n_index).sessionIndex;
% Test and Training trials are sampled from the different portation of data
nSessionData = shuffleSessionData(nDataSet(simultaneous_index), totTargets, length(totTargets));
% nSessionData = nSessionData + randn(size(nSessionData))*noise_amp/sqrt(numTrials)* addNoise(nData);
% nSessionData = normalizationDim(nSessionData, 2);
% coeffs       = coeffLDA([nSessionData; nSessionData], totTargets, totTargets);
coeffs       = coeffLDA(nSessionData, totTargets);
imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, coeffs'*coeffs);
xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
caxis([-1 1]);
axis xy;
gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
box off;
hold off;
xlabel('Time (s)')
ylabel('Time (s)')

% figure;
% 
% m = ceil(sqrt(size(nSessionData,2)));
% for nNeuron = 1:size(nSessionData,2)
%     subplot(m, m, nNeuron)
%     hold on
%     plot(squeeze(mean(nSessionData(totTargets,nNeuron,:),1)),'-b');
%     plot(squeeze(mean(nSessionData(~totTargets,nNeuron,:),1)),'-r');
%     hold off
% end


subplot(1, 3, 3)
title('Non-simultaneous Data')
hold on
random_session = randperm(length(nDataSet3D), sum(simultaneous_index));
random_index   = ones(sum(simultaneous_index),1);
for n_index    = 1:length(random_index)
    random_index(n_index) = find([nDataSet(:).sessionIndex] == ...
        nDataSet3D(random_session(n_index)).sessionIndex &  ...
        [nDataSet(:).nUnit] == nDataSet3D(random_session(n_index)).nUnit(...
        ceil(length(nDataSet3D(random_session(n_index)).nUnit)*rand())));
end
% Test and Training trials are sampled from the different portation of data
nSessionData = shuffleSessionData(nDataSet(random_index), totTargets, numTrials/2);
% nSessionData = nSessionData + randn(size(nSessionData))*noise_amp/sqrt(numTrials)* addNoise(nData);
% nSessionData = normalizationDim(nSessionData, 2);
% coeffs       = coeffLDA([nSessionData; nSessionData], totTargets, totTargets);
coeffs       = coeffLDA(nSessionData, totTargets);
imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, coeffs'*coeffs);
xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
caxis([-1 1]);
axis xy;
gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
box off;
hold off;
xlabel('Time (s)')
ylabel('Time (s)')

% figure;
% 
% m = ceil(sqrt(size(nSessionData,2)));
% for nNeuron = 1:size(nSessionData,2)
%     subplot(m, m, nNeuron)
%     hold on
%     plot(squeeze(mean(nSessionData(totTargets,nNeuron,:),1)),'-b');
%     plot(squeeze(mean(nSessionData(~totTargets,nNeuron,:),1)),'-r');
%     hold off
% end
