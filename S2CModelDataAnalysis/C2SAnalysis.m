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
    nUnitData     = nDataSet(nCell).unit_yes_trial;
    nDataSet(nCell).mcmc_yes_trial = imagingToSpike(nUnitData)/binsize;
    
    nUnitData     = nDataSet(nCell).unit_no_trial;
    nDataSet(nCell).mcmc_no_trial  = imagingToSpike(nUnitData)/binsize;
    
end

save(['S2CC2S_' DataSetList(nData).name '.mat'], 'nDataSet');