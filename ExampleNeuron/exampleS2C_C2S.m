addpath('../Func');
setDir;

load ([TempDatDir 'DataListShuffle.mat']);
nData = 1; % plot raster and psth
load([TempDatDir DataSetList(nData).name '.mat'])
spikeDataSet = nDataSet;   
load([TempDatDir 'DataListS2CModel.mat'])
params  = DataSetList(nData).params;

nDataSet = spikeDataSet;
binsize  = params.binsize;

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
    
    
    pTime         = binsize:binsize:length(paramNCells.timeSeries)*binsize;
    deconv_filter = (1 - exp(-pTime/tau_r)).* exp(-pTime/tau_d);
    C             = gallery('circul', [deconv_filter, zeros(1, length(deconv_filter))]);
    C             = C(1:length(deconv_filter), 1:length(deconv_filter));
    C             = C';
    
    nUnitData        = nCellDataSet.unit_yes_trial_linear;
    for nTrial = 1:size(nUnitData, 1)
        nUnitData(nTrial, :)        = inv(C) * nUnitData(nTrial, :)'/binsize;
    end
    nDataSet(nCell).decv_yes_trial  = nUnitData(:, 56:132);

    nUnitData        = nCellDataSet.unit_no_trial_linear;
    for nTrial = 1:size(nUnitData, 1)
        nUnitData(nTrial, :)        = inv(C) * nUnitData(nTrial, :)'/binsize;
    end
    nDataSet(nCell).decv_no_trial   = nUnitData(:, 56:132);
    
    
    nUnitData     = nCellDataSet.unit_yes_trial_linear;
    nDataSet(nCell).mcmc_yes_trial = imagingToSpike(nUnitData)/binsize;
    % nDataSet(nCell).mcmc_yes_trial = nDataSet(nCell).mcmc_yes_trial(:, 56:132);
    
    nUnitData     = nCellDataSet.unit_no_trial_linear;
    nDataSet(nCell).mcmc_no_trial  = imagingToSpike(nUnitData)/binsize;
    % nDataSet(nCell).mcmc_no_trial  = nDataSet(nCell).mcmc_no_trial(:, 56:132);
    
end

save('S2CC2S.mat', 'nDataSet');

load('S2CC2S.mat', 'nDataSet')
sim_matrix = nan(length(nDataSet), 2);

for nCell = 1:length(nDataSet)
    unit_raw  = [nDataSet(nCell).unit_yes_trial; nDataSet(nCell).unit_no_trial];
    unit_decv = [nDataSet(nCell).decv_yes_trial; nDataSet(nCell).decv_no_trial];
    unit_mcmc = [nDataSet(nCell).mcmc_yes_trial(:, 56:132); nDataSet(nCell).mcmc_no_trial(:, 56:132)];
    
    sim_matrix(nCell, 1) = corr(unit_raw(:), unit_decv(:), 'type', 'spearman');
    sim_matrix(nCell, 2) = corr(unit_raw(:), unit_mcmc(:), 'type', 'spearman');
end

figure
plot(sim_matrix(:,1), sim_matrix(:,2), 'ok', 'markersize', 2)
hold on
plot([0 1], [0 1], '--r')
xlim([0 1])
ylim([0 1])
xlabel('Similarity ephys pattern')
ylabel('Similarity ephys pattern')
setPrint(8, 6, 'S2C_C2S2')

