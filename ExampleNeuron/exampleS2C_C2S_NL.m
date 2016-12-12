addpath('../Func');
setDir;

load ([TempDatDir 'DataListShuffle.mat']);
nData = 1; % plot raster and psth
load([TempDatDir DataSetList(nData).name '.mat'])
spikeDataSet = nDataSet;
load([TempDatDir 'DataListS2CModel.mat'])
nData = 2;
params  = DataSetList(nData).params;
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');
g = @(p,x) p(1) + p(2)./ (1 + exp((p(3)-x)*p(4)));
inv_g = @(p, x) p(3) - 1/p(4) * log(p(2)./(x-p(1)) - 1);

int_noise = [2.0, 2.0]; % 6s-AAV, GP4.3
ext_noise = [0.15, 0.45];

nDataSet = spikeDataSet;
binsize  = params.binsize;
nData = 2;

for nCell = 1:length(nDataSet)
    tau_d = params.tau_d(nCell);
    tau_r = params.tau_r(nCell);
    paramNCells = params;
    paramNCells.Fm = 1;
    paramNCells.K  = 1;
    paramNCells.n  = 1;
    paramNCells.tau_r = tau_r;
    paramNCells.tau_d = tau_d;
    paramNCells.intNoise = 0;
    paramNCells.extNoise = 0;
    paramNCells.timeSeries = -100*binsize:binsize:5;
    nCellDataSet  = getFakeCaImagingData(spikeDataSet(nCell), paramNCells);
    % I only put external noise here
    yesNoise = randn(size(nCellDataSet.unit_yes_trial))*ext_noise(nData);
    noNoise  = randn(size(nCellDataSet.unit_no_trial))*ext_noise(nData);
    yesIntNoise = randn(size(nCellDataSet.unit_yes_trial))*int_noise(nData);
    noIntNoise  = randn(size(nCellDataSet.unit_no_trial))*int_noise(nData);
    
    % a parameter from random cell
    randCell      = ceil(rand*length(nDataSet));
% %     tau_d         = params.tau_d(randCell);
% %     tau_r         = params.tau_r(randCell);
    tau_d         = params.tau_d(nCell);
    tau_r         = params.tau_r(nCell);
    pTime         = binsize:binsize:length(paramNCells.timeSeries)*binsize;
    deconv_filter = (1 - exp(-pTime/tau_r)).* exp(-pTime/tau_d);
    C             = gallery('circul', [deconv_filter, zeros(1, length(deconv_filter))]);
    C             = C(1:length(deconv_filter), 1:length(deconv_filter));
    C             = C';
    invC          = inv(C);
    
    nUnitData     = nCellDataSet.unit_yes_trial_linear;
    nUnitData     = nUnitData + yesIntNoise;
    nUnitData(nUnitData<0) = 0;
    param         = squeeze(nlParams(2, nCell, :));
    paramInv      = squeeze(nlParams(2, nCell, :));
% % %     paramInv      = squeeze(nlParams(2, randCell, :));
% % %     paramInv = squeeze(nlParams(2, nCell, :));
    nUnitData     = g(param, nUnitData) + yesNoise;
    yesUnitData   = nUnitData;
    
    nUnitData     = nCellDataSet.unit_no_trial_linear;
    nUnitData     = nUnitData + noIntNoise;
    nUnitData(nUnitData<0) = 0;
    nUnitData     = g(param, nUnitData) + noNoise;
    noUnitData    = nUnitData;
% %     
% %     minData       = min([mean(yesUnitData,1), mean(noUnitData,1)]);
% %     maxData       = max([mean(yesUnitData,1), mean(noUnitData,1)]);
% %     
% %     paramInv(1)   = minData; % param(1);
% %     paramInv(2)   = maxData; % param(2);
    
    % average across trial
    nUnitData = yesUnitData;
%     nUnitData = mean(yesUnitData, 1);
    % cutting the top and bottom 5% signal
    per_cent = 0.02;
    nUnitData(nUnitData > paramInv(1) + paramInv(2)*(1-per_cent)) = paramInv(1) + paramInv(2)*(1-per_cent); 
    nUnitData(nUnitData < paramInv(1) + paramInv(2)*per_cent) = paramInv(1) + paramInv(2)*per_cent; 
    nUnitData = inv_g(paramInv, nUnitData);
    
    for nTrial = 1:size(nUnitData, 1)
%         XN            = nUnitData(nTrial, :);
%         xdMODWT       = wden(XN,'modwtsqtwolog','s','mln',4,'sym4');
        xdMODWT         = nUnitData(nTrial, :);
        nUnitData(nTrial, :)        = invC * xdMODWT'/binsize;
    end
    nDataSet(nCell).decv_yes_trial  = nUnitData(:, 56:132);


    % average across trial
    nUnitData = noUnitData;
%     nUnitData = mean(noUnitData, 1);
    % cutting the top and bottom 5% signal
    nUnitData(nUnitData > paramInv(1) + paramInv(2)*(1-per_cent)) = paramInv(1) + paramInv(2)*(1-per_cent); 
    nUnitData(nUnitData < paramInv(1) + paramInv(2)*per_cent) = paramInv(1) + paramInv(2)*per_cent; 
    nUnitData = inv_g(paramInv, nUnitData);    
    
    for nTrial = 1:size(nUnitData, 1)
%         XN            = nUnitData(nTrial, :);
%         xdMODWT       = wden(XN,'modwtsqtwolog','s','mln',4,'sym4');
        xdMODWT         = nUnitData(nTrial, :);
        nUnitData(nTrial, :)        = invC * xdMODWT'/binsize;
    end
    nDataSet(nCell).decv_no_trial   = nUnitData(:, 56:132);
    
end

save('DeconvGPResults/S2CC2SNLPreciseSingleTrial.mat', 'nDataSet');

% % % load('DeconvGPResults/S2CC2SNLRandomNLSingleTrial.mat', 'nDataSet')
% % % sim_matrix = nan(length(nDataSet), 1);
% % % 
% % % ephyAct  = nan(length(nDataSet), 1);
% % % 
% % % for nCell = 1:length(nDataSet)
% % %     unit_raw  = [mean(nDataSet(nCell).unit_yes_trial, 1); mean(nDataSet(nCell).unit_no_trial, 1)];
% % %     unit_decv = [mean(nDataSet(nCell).decv_yes_trial, 1); mean(nDataSet(nCell).decv_no_trial, 1)];
% % %     sim_matrix(nCell) = corr(unit_raw(:), unit_decv(:), 'type', 'spearman');
% % %     ephyAct(nCell) = mean(unit_raw(:));
% % % end
% % % 
% % % load('DeconvGPResults/S2CC2SMCMCSingleTrial_FineTunedModeled_6s_AAV.mat', 'nDataSet')
% % % sim_matrix = nan(length(nDataSet), 1);
% % % for nCell = 1:length(nDataSet)
% % %     unit_raw  = [mean(nDataSet(nCell).unit_yes_trial, 1); mean(nDataSet(nCell).unit_no_trial, 1)];
% % %     unit_decv = [mean(nDataSet(nCell).mcmc_yes_trial, 1); mean(nDataSet(nCell).mcmc_no_trial, 1)];
% % %     sim_matrix(nCell) = corr(unit_raw(:), unit_decv(:), 'type', 'spearman');
% % % end
% % % 
% % % 
% % % load('DeconvResults/S2CC2SNLPreciseAve.mat', 'nDataSet')
% % % sim_matrixRandom = nan(length(nDataSet), 1);
% % % for nCell = 1:length(nDataSet)
% % %     unit_raw  = [mean(nDataSet(nCell).unit_yes_trial, 1); mean(nDataSet(nCell).unit_no_trial, 1)];
% % %     unit_decv = [mean(nDataSet(nCell).decv_yes_trial, 1); mean(nDataSet(nCell).decv_no_trial, 1)];
% % %     sim_matrixRandom(nCell) = corr(unit_raw(:), unit_decv(:), 'type', 'spearman');
% % % end
% % % 
% % % 
% % % figure
% % % plot(sim_matrixRandom, sim_matrix, 'ok')
% % % hold on
% % % plot([-0.2, 1], [-0.2, 1], '--r')
% % % xlabel('Similarity index -- random deconv + single trial')
% % % ylabel('Similarity index -- precise deconv + average trial')
% % % xlim([-0.2 1.01])
% % % ylim([-0.2 1.01])
% % % setPrint(8, 6, 'Similarity_index')
% % % 
% % % [p, h] = ttest(sim_matrixRandom, sim_matrix)