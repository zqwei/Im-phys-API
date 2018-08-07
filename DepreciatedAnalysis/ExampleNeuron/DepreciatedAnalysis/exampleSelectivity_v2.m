addpath('../Func');
setDir;

load ([TempDatDir 'DataListShuffle.mat']);
nData = 1; % plot raster and psth
load([TempDatDir DataSetList(nData).name '_old.mat'])
numT       = size(nDataSet(1).unit_yes_trial,2);
params     = DataSetList(1).params;

sigma                         = 0.15 / params.binsize; % 200 ms
filterLength                  = 11;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse); 


actMat     = nan(length(nDataSet),numT*2);
for nUnit = 1: length(nDataSet)
    unitYesAct                = mean(getGaussianPSTH (filterInUse, nDataSet(nUnit).unit_yes_trial, 2));
    unitNoAct                 = mean(getGaussianPSTH (filterInUse, nDataSet(nUnit).unit_no_trial, 2));
    actMat(nUnit, 1:numT)     = unitYesAct;
    actMat(nUnit, numT+1:end) = unitNoAct;
end


load ([TempDatDir 'DataListS2CModel.mat']);
nData = 4; % plot s2c model
load([TempDatDir DataSetList(nData).name '.mat'])
s2cActMat    = nan(length(nDataSet),numT*2);
for nUnit = 1: length(nDataSet)
    unitYesAct                = mean(getGaussianPSTH (filterInUse, nDataSet(nUnit).unit_yes_trial, 2));
    unitNoAct                 = mean(getGaussianPSTH (filterInUse, nDataSet(nUnit).unit_no_trial, 2));
    s2cActMat(nUnit, 1:numT)     = unitYesAct;
    s2cActMat(nUnit, numT+1:end) = unitNoAct;
end

corrAct = corr(actMat');
corrSAct = corr(s2cActMat');
corrAct = abs(corrAct);
corrSAct = abs(corrSAct);
plot(corrSAct - corrAct)