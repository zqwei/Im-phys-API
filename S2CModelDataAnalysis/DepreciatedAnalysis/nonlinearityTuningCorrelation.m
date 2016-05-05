addpath('../Func');
setDir;

load ([TempDatDir 'DataListShuffle.mat']);
nData = 4;
load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
numTrials        = 100;
firingRates      = generateDPCAData(nDataSet, numTrials);
caFiringRates    = nanmean(firingRates, ndims(firingRates));
caFiringRates    = [squeeze(caFiringRates(:, 1, :)), squeeze(caFiringRates(:, 2, :))];


nData = 1;
load([TempDatDir DataSetList(nData).name '.mat'])
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
spikeFiringRates       = nanmean(firingRates, ndims(firingRates));
spikeFiringRates       = [squeeze(spikeFiringRates(:, 1, :)), squeeze(spikeFiringRates(:, 2, :))];
% find the spearman correlation
rhoMat                 = corr(caFiringRates', spikeFiringRates', 'type', 'Spearman');
[maxValue, maxIndex]   = max(rhoMat, [], 1);

mapFiringRates         = caFiringRates(maxIndex, :);

plot()