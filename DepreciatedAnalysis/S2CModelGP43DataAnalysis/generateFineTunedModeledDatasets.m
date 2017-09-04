%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating new data using fine tuned params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% one should run:
% 1.
% generateModeledDatasets
% 2.
% nonlinearityTuning
% 
% to first obtain the linear filtered signal
% to second obtain nonlinearity and noise

addpath('../Func');
setDir;
load('FineTunedNLParams.mat');
load ([TempDatDir 'DataListS2CGP43Model.mat']);

% nonlinear function
% g = @(p,x) p(1) + p(2)./ (1 + (p(3)./x).^p(4));
g = @(p,x) p(1) + p(2)./ (1 + 10.^((p(3)-x)*p(4)));

numTrial = 50;

ext_noise = 0.90;
% numDataList = length(DataSetList);
numDataList = 1;

for nData = 1:numDataList
    
    DataSetList(numDataList+nData) = DataSetList(nData);
    DataSetList(numDataList+nData).name = ['FineTuned' DataSetList(nData).name];
    
    load([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');
    for nUnit  = 1:length(nDataSet)
        param  = nlParams(nUnit, :);
        yesNoise = randn(numTrial, 77)*ext_noise; %squeeze(noiseRates{nData}(nUnit, 1, :, :))';
        noNoise  = randn(numTrial, 77)*ext_noise; %squeeze(noiseRates{nData}(nUnit, 2, :, :))';
        nDataSet(nUnit).unit_yes_trial = g(param, nDataSet(nUnit).unit_yes_trial_linear) + ...
                                            yesNoise(randpermLargeK(numTrial, size(nDataSet(nUnit).unit_yes_trial, 1)), :);
        nDataSet(nUnit).unit_no_trial  = g(param, nDataSet(nUnit).unit_no_trial_linear) + ...
                                            noNoise(randpermLargeK(numTrial, size(nDataSet(nUnit).unit_no_trial, 1)), :);
        nDataSet(nUnit).unit_yes_error = g(param, nDataSet(nUnit).unit_yes_error_linear) + ...
                                            yesNoise(randpermLargeK(numTrial, size(nDataSet(nUnit).unit_yes_error, 1)), :);
        nDataSet(nUnit).unit_no_error  = g(param, nDataSet(nUnit).unit_no_error_linear) + ...
                                            noNoise(randpermLargeK(numTrial, size(nDataSet(nUnit).unit_no_error, 1)), :);
    end
    save([TempDatDir DataSetList(numDataList+nData).name '.mat'], 'nDataSet');     
end

save([TempDatDir 'DataListS2CGP43Model.mat'], 'DataSetList');