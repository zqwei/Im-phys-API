addpath('../Func');
setDir;
load([TempDatDir 'FineTuned6fNLParams.mat'], 'nlParams');
load([TempDatDir 'DataListS2C6fModel.mat']);

% nonlinear function
g = @(p,x) p(1) + p(2)./ (1 + exp((p(3)-x)*p(4)));

int_noise = 0.95;
ext_noise = 0.40;

nData = 1;
medianParam = nanmedian(nlParams);
load([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');
nlParams_   = nlParams;
for nUnit   = 1:length(nDataSet)
    param   = squeeze(nlParams(nUnit, :));
    if sum(isnan(param)) > 0
        param = medianParam;
        nlParams_(nUnit, :) = medianParam;
    end
    yesNoise = randn(size(nDataSet(nUnit).unit_yes_trial))*ext_noise(nData);
    noNoise  = randn(size(nDataSet(nUnit).unit_no_trial))*ext_noise(nData);
    yesIntNoise = randn(size(nDataSet(nUnit).unit_yes_trial))*int_noise(nData);
    noIntNoise  = randn(size(nDataSet(nUnit).unit_no_trial))*int_noise(nData);
    unit_yes_trial_linear = nDataSet(nUnit).unit_yes_trial_linear + yesIntNoise;
    unit_no_trial_linear  = nDataSet(nUnit).unit_no_trial_linear + noIntNoise;
    unit_yes_trial_linear(unit_yes_trial_linear<0) = 0;
    unit_no_trial_linear(unit_no_trial_linear<0)   = 0;
    nDataSet(nUnit).unit_yes_trial = g(param, unit_yes_trial_linear) + yesIntNoise;
    nDataSet(nUnit).unit_no_trial  = g(param, unit_no_trial_linear) + noIntNoise;
end

caDataSet              = nDataSet;
tau_r                  = DataSetList(nData).params.tau_r;
tau_d                  = DataSetList(nData).params.tau_d;
params                 = DataSetList(nData).params;
nDataSet               = getFakeSpikeNLDeconvData(caDataSet, tau_r, tau_d, nlParams, params);
% figure
% plotMeanActivityImagescRasterOnlyPositivePeak(nDataSet, DataSetList(1).params, [], []); 
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
fName    = 'DeconvNL_Ca_Fast_SShort_Delay';
save(['S2CC2S_' fName '.mat'], 'nDataSet');

truncatedNormal        = truncate(makedist('Normal'), -0.9, 1.5);
params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
std_r                  = 0.0222;
median_r               = 0.0213;
std_d                  = 0.5390;
median_d               = 0.5898;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;
nDataSet               = getFakeSpikeDeconvData(caDataSet, params.tau_r, params.tau_d, params);  
% figure
% plotMeanActivityImagescRasterOnlyPositivePeak(nDataSet, DataSetList(1).params, [], []); 
fName    = 'Deconv_Ca_Fast_SShort_Delay';
save(['S2CC2S_' fName '.mat'], 'nDataSet');