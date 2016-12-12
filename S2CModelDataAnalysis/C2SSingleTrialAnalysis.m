addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CModel.mat']);

nData = 1; % plot raster and psth
load([TempDatDir DataSetList(nData).name '.mat'])
params  = DataSetList(nData).params;
binsize  = params.binsize;


nData = 3;
load([TempDatDir DataSetList(nData).name '.mat'])
for nCell = 1:length(nDataSet)
    
    numYesTrial      = size(nDataSet(nData).unit_yes_trial, 1);
    numNoTrial       = size(nDataSet(nData).unit_no_trial, 1);
    numT             = size(nDataSet(nData).unit_yes_trial, 2);
    
    nUnitData        = [nDataSet(nCell).unit_yes_trial; nDataSet(nCell).unit_no_trial];
    nUnitData        = nUnitData';
    
    fastData         = imagingToSpike(nUnitData(:)');
    fastData         = reshape(fastData, numT, numYesTrial + numNoTrial);
    
    nDataSet(nCell).mcmc_yes_trial = fastData(:, 1:numYesTrial)'/binsize;    
    nDataSet(nCell).mcmc_no_trial  = fastData(:, 1+numYesTrial:end)'/binsize;
    
end

save(['S2CC2SMCMCSingleTrial_' DataSetList(nData).name '.mat'], 'nDataSet');