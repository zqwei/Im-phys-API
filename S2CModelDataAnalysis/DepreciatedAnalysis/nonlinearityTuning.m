addpath('../Func');
setDir;

load ([TempDatDir 'DataListShuffle.mat']);
% nData = 4;
% load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
% params      = DataSetList(nData).params;
% caContraIndex = false(length(nDataSet), 1);
% yesActMat   = nan(length(nDataSet), length(params.timeSeries));
% noActMat    = nan(length(nDataSet), length(params.timeSeries));
% timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
% for nUnit   = 1:length(nDataSet)
%     yesTrial = mean(nDataSet(nUnit).unit_yes_trial);
%     noTrial  = mean(nDataSet(nUnit).unit_no_trial);
%     yesActMat(nUnit, :)  = yesTrial;
%     noActMat(nUnit, :)   = noTrial;
%     caContraIndex(nUnit)   = sum(noTrial(timePoints(2):end))>sum(yesTrial(timePoints(2):end));
% end
% 
% figure;
% hold on
% shadedErrorBar(params.timeSeries, mean(noActMat(caContraIndex,:)), sem(noActMat(caContraIndex,:)),'-b')
% shadedErrorBar(params.timeSeries, mean(yesActMat(caContraIndex,:)), sem(yesActMat(caContraIndex,:)),'-r')
% title('contra neuron only')
% 
% figure;
% hold on
% shadedErrorBar(params.timeSeries, mean(noActMat(~caContraIndex,:)), sem(noActMat(~caContraIndex,:)),'-b')
% shadedErrorBar(params.timeSeries, mean(yesActMat(~caContraIndex,:)), sem(yesActMat(~caContraIndex,:)),'-r')
% title('ipsi neuron only')
% 
% refcontraMat = [mean(yesActMat(caContraIndex,:)), mean(noActMat(caContraIndex,:))];
% refipsiMat   = [mean(yesActMat(~caContraIndex,:)), mean(noActMat(~caContraIndex,:))];
% save('refMat', 'refcontraMat', 'refipsiMat')

nData = 1;
load([TempDatDir DataSetList(nData).name '.mat'])
% params      = DataSetList(nData).params;
% spikeContraIndex = false(length(nDataSet), 1);
% cellType    = [DataSetList(nData).cellinfo.cellType]';
% yesActMat   = nan(length(nDataSet), length(params.timeSeries));
% noActMat    = nan(length(nDataSet), length(params.timeSeries));
% timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
% for nUnit   = 1:length(nDataSet)
%     yesTrial = mean(nDataSet(nUnit).unit_yes_trial);
%     noTrial  = mean(nDataSet(nUnit).unit_no_trial);
%     yesActMat(nUnit, :)  = yesTrial;
%     noActMat(nUnit, :)   = noTrial;
%     spikeContraIndex(nUnit)   = sum(noTrial(timePoints(2):end))>sum(yesTrial(timePoints(2):end));
% end
% save('spikeContraIndex', 'spikeContraIndex')
% valid_depth = [nDataSet.depth_in_um];
% valid_depth = valid_depth < 750 & valid_depth > 140;
% valid_depth = valid_depth';

load('refMat')
load('spikeContraIndex')
load('ParamsFitCells_S2CModel_Sim.mat');
S2Cparams   = params;
clear params;
idTauD      = 1;
idN         = 1;
idK         = 1;
params      = DataSetList(nData).params;
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
params.Fm              = ones(length(nDataSet), 1)*21.3240;
params.K               = ones(length(nDataSet), 1)*S2Cparams(1).K;
params.n               = ones(length(nDataSet), 1)*S2Cparams(1).n;
params.tau_r           = ones(length(nDataSet), 1)*S2Cparams(1).tau_r;
params.tau_d           = ones(length(nDataSet), 1)*S2Cparams(idTauD).tau_d;
params.Fm              = params.Fm;
params.tau_r           = params.tau_r;
params.tau_d           = params.tau_d*1.0;
params.intNoise        = 1.5;
params.extNoise        = 0;
oldDataSet             = getFakeCaImagingData(nDataSet, params);
nDataSet               = generateNonLinearDataFromLinear(oldDataSet, params);

numTrials              = 100;
firingRates            = generateDPCAData(nDataSet, numTrials);
firingRatesAverage     = nanmean(firingRates, ndims(firingRates));
firingRatesAverage     = [squeeze(firingRatesAverage(:, 1, :)), squeeze(firingRatesAverage(:, 2, :));];
linearcontraMat        = mean(firingRatesAverage(spikeContraIndex & valid_depth, :));
linearipsiMat          = mean(firingRatesAverage(~spikeContraIndex & valid_depth, :));
[param,~]=sigm_fit([linearcontraMat'; linearipsiMat'], [refcontraMat'; refipsiMat'], [],[-0.07, 0.24, 1, 1.5*log(10)]);
param(2) = param(2) - param(1);
f = @(p,x) p(1) + p(2)./ (1 + 10.^((p(3)-x)*p(4)));

% 
% firingRatesYes         = squeeze(firingRates(:,1,:,:));
% firingRatesNo          = squeeze(firingRates(:,2,:,:));
% firingRatesAll         = [permute(firingRatesYes, [2 1 3]); permute(firingRatesNo, [2 1 3])];
% firingRatesIpsi        = firingRatesAll(:, ~spikeContraIndex, :);
% firingRatesContra      = firingRatesAll(:, spikeContraIndex, :);

fp1 = @(p, x) ones(size(x));
% fp2 = @(p, x) 1./(10.^(p(4)*(p(3)-x)) + 1);
% fp3 = @(p, x) -(10.^(p(4)*(p(3)-x))*p(2)*p(4)*log(10))./((10.^(p(4)*(p(3)-x)) + 1).^2);
% fp4 = @(p, x) 10.^(p(4)*(p(3)-x))*p(2)*log(10).*(x - p(3))./((10.^(p(4)*(p(3)-x)) + 1).^2);
% 
% ita0 = 1;
% 
% for nError             = 1:10000
%     nlfiringRates      = f(param, firingRates);
%     firingRatesAverage = nanmean(nlfiringRates, ndims(nlfiringRates));
%     yesActMat          = squeeze(firingRatesAverage(:, 1, :));
%     noActMat           = squeeze(firingRatesAverage(:, 2, :));
%     nlcontraMat        = [mean(yesActMat(spikeContraIndex,:)), mean(noActMat(spikeContraIndex,:))];
%     nlipsiMat          = [mean(yesActMat(~spikeContraIndex,:)), mean(noActMat(~spikeContraIndex,:))];
%     errorMat           = [nlcontraMat'; nlipsiMat']- [refcontraMat'; refipsiMat'];
%     errorSqr           = mean(errorMat.^2);
%     disp(errorSqr)
%     fp1Mat             = fp1(param, firingRates);
%     fp2Mat             = fp2(param, firingRates);
%     fp3Mat             = fp3(param, firingRates);
%     fp4Mat             = fp4(param, firingRates);
% 
%     ita                = ita0/(20+nError)*20;    
%     
%     meanFP1            = nanmean(fp1Mat, ndims(fp1Mat));
%     yesActMat          = squeeze(meanFP1(:, 1, :));
%     noActMat           = squeeze(meanFP1(:, 2, :));
%     nlcontraMat        = [mean(yesActMat(spikeContraIndex,:)), mean(noActMat(spikeContraIndex,:))];
%     nlipsiMat          = [mean(yesActMat(~spikeContraIndex,:)), mean(noActMat(~spikeContraIndex,:))];
%     updateFP1          = mean(errorMat.*([nlcontraMat'; nlipsiMat']));
%     param(1)           = param(1) - updateFP1 * ita;
% 
%     meanFP2            = nanmean(fp2Mat, ndims(fp2Mat));
%     yesActMat          = squeeze(meanFP2(:, 1, :));
%     noActMat           = squeeze(meanFP2(:, 2, :));
%     nlcontraMat        = [mean(yesActMat(spikeContraIndex,:)), mean(noActMat(spikeContraIndex,:))];
%     nlipsiMat          = [mean(yesActMat(~spikeContraIndex,:)), mean(noActMat(~spikeContraIndex,:))];
%     updateFP2          = mean(errorMat.*([nlcontraMat'; nlipsiMat']));
%     param(2)           = param(2) - updateFP2 * ita;
% 
%     meanFP3            = nanmean(fp3Mat, ndims(fp3Mat));
%     yesActMat          = squeeze(meanFP3(:, 1, :));
%     noActMat           = squeeze(meanFP3(:, 2, :));
%     nlcontraMat        = [mean(yesActMat(spikeContraIndex,:)), mean(noActMat(spikeContraIndex,:))];
%     nlipsiMat          = [mean(yesActMat(~spikeContraIndex,:)), mean(noActMat(~spikeContraIndex,:))];
%     updateFP3          = mean(errorMat.*([nlcontraMat'; nlipsiMat']));
%     param(3)           = param(3) - updateFP3 * ita;
%     
% 
%     meanFP4            = nanmean(fp4Mat, ndims(fp4Mat));
%     yesActMat          = squeeze(meanFP4(:, 1, :));
%     noActMat           = squeeze(meanFP4(:, 2, :));
%     nlcontraMat        = [mean(yesActMat(spikeContraIndex,:)), mean(noActMat(spikeContraIndex,:))];
%     nlipsiMat          = [mean(yesActMat(~spikeContraIndex,:)), mean(noActMat(~spikeContraIndex,:))];
%     updateFP4          = mean(errorMat.*([nlcontraMat'; nlipsiMat']));
%     param(4)           = param(4) - updateFP4 * ita;
%     
% end


nlfiringRates = f(param, firingRates);
firingRatesAverage     = nanmean(nlfiringRates, ndims(nlfiringRates));
% firingRatesAverage     = nanmean(firingRates, ndims(firingRates));
% firingRatesAverage     = f(param, firingRatesAverage);
yesActMat              = squeeze(firingRatesAverage(:, 1, :));
noActMat               = squeeze(firingRatesAverage(:, 2, :));

figure;
hold on
shadedErrorBar(params.timeSeries, mean(noActMat(spikeContraIndex & valid_depth,:)), sem(noActMat(spikeContraIndex & valid_depth,:)),'-b')
shadedErrorBar(params.timeSeries, mean(yesActMat(spikeContraIndex & valid_depth,:)), sem(yesActMat(spikeContraIndex & valid_depth,:)),'-r')
title('contra neuron only')

figure;
hold on
shadedErrorBar(params.timeSeries, mean(noActMat(~spikeContraIndex & valid_depth,:)), sem(noActMat(~spikeContraIndex & valid_depth,:)),'-b')
shadedErrorBar(params.timeSeries, mean(yesActMat(~spikeContraIndex & valid_depth,:)), sem(yesActMat(~spikeContraIndex & valid_depth,:)),'-r')
title('ipsi neuron only')

% figure;
% hold on
% plot(f(param, linearcontraMat(1:77)'), '--b')
% plot(refcontraMat(1:77)','-b')
% plot(f(param, linearcontraMat(78:end)'), '--r')
% plot(refcontraMat(78:end)','-r')
% figure;
% hold on
% plot(f(param, linearipsiMat(1:77)'), '--b')
% plot(refipsiMat(1:77)','-b')
% plot(f(param, linearipsiMat(78:end)'), '--r')
% plot(refipsiMat(78:end)','-r')

