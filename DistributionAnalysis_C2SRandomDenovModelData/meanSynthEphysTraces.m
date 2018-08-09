%
% Compute distribution of selectivity for each single neuron
% 
% Here we drop the test of parameter of tau_r
% 
% 
% -------------------------------------------------------------------------
% version 1.0
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
numFold = 30;

% short Ca Fast GP517
nData  = 2;
params = DataSetList(nData).params;
std_r                  = 0.0222;
median_r               = 0.0212;
std_d                  = 0.3447;
median_d               = 0.5656;
load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
per_list               = 0.02:0.01:0.98;
tau_d_list             = icdf('Normal', per_list, 0, 1) * std_d + median_d;

for nTau              = 1:length(tau_d_list)
    if tau_d_list(nTau) > 0
        spikeDataSet      = getSyntheticSpikeDeconvDataSimpleVersion(nDataSet, median_r, tau_d_list(nTau), params);   %#ok<*NASGU>
        save([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'spikeDataSet');
        clear spikeDataSet
    end
end


% short Ca slow
nData  = 3;
params = DataSetList(nData).params;
std_r                  = 0.0375;
median_r               = 0.0927;
std_d                  = 0.5374;
median_d               = 1.2294;
load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
per_list               = 0.02:0.01:0.98;
tau_d_list             = icdf('Normal', per_list, 0, 1) * std_d + median_d;

for nTau              = 1:length(tau_d_list)
    if tau_d_list(nTau) > 0
        spikeDataSet      = getSyntheticSpikeDeconvDataSimpleVersion(nDataSet, median_r, tau_d_list(nTau), params);   %#ok<*NASGU>
        save([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'spikeDataSet');
        clear spikeDataSet
    end
end


% short Ca slow virus
nData  = 4;
params = DataSetList(nData).params;
std_r                  = 0.0246;
median_r               = 0.0505;
std_d                  = 0.4588;
median_d               = 1.7064;
load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
per_list               = 0.02:0.01:0.98;
tau_d_list             = icdf('Normal', per_list, 0, 1) * std_d + median_d;

for nTau              = 1:length(tau_d_list)
    if tau_d_list(nTau) > 0
        spikeDataSet      = getSyntheticSpikeDeconvDataSimpleVersion(nDataSet, median_r, tau_d_list(nTau), params);   %#ok<*NASGU>
        save([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'spikeDataSet');
        clear spikeDataSet
    end
end


% long Ca fast
nData  = 5;
params = DataSetList(nData).params;
std_r                  = 0.0222;
median_r               = 0.0212;
std_d                  = 0.3447;
median_d               = 0.5656;
load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
per_list               = 0.02:0.01:0.98;
tau_d_list             = icdf('Normal', per_list, 0, 1) * std_d + median_d;

for nTau              = 1:length(tau_d_list)
    if tau_d_list(nTau) > 0
        spikeDataSet      = getSyntheticSpikeDeconvDataSimpleVersion(nDataSet, median_r, tau_d_list(nTau), params);   %#ok<*NASGU>
        save([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'spikeDataSet');
        clear spikeDataSet
    end
end


% long Ca slow
nData  = 6;
params = DataSetList(nData).params;
std_r                  = 0.0375;
median_r               = 0.0927;
std_d                  = 0.5374;
median_d               = 1.2294;
load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
per_list               = 0.02:0.01:0.98;
tau_d_list             = icdf('Normal', per_list, 0, 1) * std_d + median_d;

for nTau              = 1:length(tau_d_list)
    if tau_d_list(nTau) > 0
        spikeDataSet      = getSyntheticSpikeDeconvDataSimpleVersion(nDataSet, median_r, tau_d_list(nTau), params);   %#ok<*NASGU>
        save([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'spikeDataSet');
        clear spikeDataSet
    end
end

% 6f data
nData  = 10;
params = DataSetList(nData).params;
std_r                  = 0.0222;
median_r               = 0.0212;
std_d                  = 0.3447;
median_d               = 0.5656;
load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
per_list               = 0.02:0.01:0.98;
tau_d_list             = icdf('Normal', per_list, 0, 1) * std_d + median_d;
for nTau              = 1:length(tau_d_list)
    if tau_d_list(nTau) > 0
        spikeDataSet      = getSyntheticSpikeDeconvDataSimpleVersion(nDataSet, median_r, tau_d_list(nTau), params);   %#ok<*NASGU>
        save([TempDatDir 'directDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'spikeDataSet');
        clear spikeDataSet
    end
end