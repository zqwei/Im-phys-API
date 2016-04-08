%
% A code to check the influence factors among (baseline firing rate, rising
% time, decaying time, internal noise, n and k) to the generated calcium
% activity profiles.
%
% for this code, only activity profile for example neurons are generated.
% 
% for the easy of illustration in model case, I only show a change of
% parameter by two levels
%
% a more detailed comparison is based on those of 10 levels of changes
% 
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('../Func/')
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
load ([TempDatDir 'ParamsFitCells_S2CModel_Sim.mat'], 'params');
paramsSim = params;
clear params;
nData = 1;
load([TempDatDir DataSetList(nData).name '.mat'])

Fm              = 21.3240;
intNoise        = 1.5;
extNoise        = 0;
nGroup          = 2;
nNeuron         = 30;
spikeDataSet           = nDataSet(nNeuron);
params                 = DataSetList(nData).params;
params.Fm              = Fm;
params.K               = paramsSim(1, nGroup).K;
params.n               = paramsSim(1, nGroup).n;
params.tau_r           = paramsSim(1, nGroup).tau_r;
params.tau_d           = paramsSim(1, nGroup).tau_d;
params.intNoise        = intNoise;
params.extNoise        = extNoise;   


nUnitData       = spikeDataSet.unit_yes_trial;
yesTrialRate    = mean(nUnitData, 1);
nUnitData       = spikeDataSet.unit_no_trial;
noTrialRate     = mean(nUnitData, 1);
numYesTrial     = size(nUnitData, 1);
numNoTrial      = size(nUnitData, 1);
timeBins = params.timeSeries;
binSize  = params.binsize;

nFactor         = 5;

[spkTime, spkCounts] = NHpoisson(yesTrialRate * nFactor, timeBins, binSize, numYesTrial);
spikeDataSet.unit_yes_trial          = spkCounts;
spikeDataSet.unit_yes_trial_spk_time = spkTime;
[spkTime, spkCounts] = NHpoisson(noTrialRate * nFactor, timeBins, binSize, numNoTrial);
spikeDataSet.unit_no_trial           = spkCounts;
spikeDataSet.unit_no_trial_spk_time  = spkTime;


figure;

subplot(2, 2, 1)
spkTimes{1}    = spikeDataSet.unit_yes_trial_spk_time;
spkTimes{2}    = spikeDataSet.unit_no_trial_spk_time;
oldBlue         = [     0         0       0.7];
oldRed          = [    0.7        0         0];
color_index    = [oldRed; oldBlue];    
spkLoc         = 0;
for nPlot            = 1:2
    hold on;
    for nTrial       = 1:length(spkTimes{nPlot})
        spikeTimes   = spkTimes{nPlot}{nTrial};
        spikeTrial   = ones(length(spikeTimes), 1) * (nTrial + spkLoc);
        plot(spikeTimes, spikeTrial, '.', 'color', color_index(nPlot, :));
    end
    spkLoc     = spkLoc + length(spkTimes{nPlot}) + 3;
end
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0);
hold off;
ylim([1 spkLoc])
xlim([params.timeSeries(1) params.timeSeries(end)]);
axis off


subplot(2, 2, 3)
oldBlue         = [     0         0       0.7];
oldRed          = [    0.7        0         0];
color_index     = [oldRed; oldBlue];
hold on;
nUnitData       = spikeDataSet.unit_yes_trial;
yesTrialRate    = mean(nUnitData, 1);
stdYesTrial     = std(nUnitData, 1)/sqrt(size(nUnitData, 1));
nUnitData       = spikeDataSet.unit_no_trial;
noTrialRate     = mean(nUnitData, 1);
stdNoTrial      = std(nUnitData, 1)/sqrt(size(nUnitData, 1));
shadedErrorBar(params.timeSeries, yesTrialRate, stdYesTrial, ...
    {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
shadedErrorBar(params.timeSeries, noTrialRate, stdNoTrial, ...
    {'k-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
ylabel('Firing rate (Hz)');
xlabel('Time (s)');
xlim([params.timeSeries(1) params.timeSeries(end)]);
set(gca, 'TickDir', 'out')    



subplot(2, 2, 2)

nDataSet             = getFakeCaImagingData(spikeDataSet, params);
nYesData             = nDataSet.unit_yes_trial;
nNoData              = nDataSet.unit_no_trial;
stdYesTrial          = std(nYesData, 1)/sqrt(size(nYesData, 1));
stdNoTrial           = std(nNoData, 1)/sqrt(size(nNoData, 1));
nYesData             = mean(nYesData, 1);
nNoData              = mean(nNoData, 1);
maxDff               = max([nYesData, nNoData]);
minDff               = min([nYesData, nNoData]);
diffDff              = maxDff - minDff;        
nYesData             = (nYesData - minDff)/diffDff;
nNoData              = (nNoData - minDff)/diffDff;
stdYesTrial          = stdYesTrial/diffDff;
stdNoTrial           = stdNoTrial/diffDff;
hold on
shadedErrorBar(params.timeSeries, nYesData, stdYesTrial, ...
    {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
shadedErrorBar(params.timeSeries, nNoData, stdNoTrial, ...
    {'k-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
ylabel('Normalized dff (Hz)');
xlabel('Time (s)');
ylim([0 1])
xlim([params.timeSeries(1) params.timeSeries(end)]);
set(gca, 'TickDir', 'out')    
title('S2C model')


subplot(2, 2, 4)
nDataSet             = getFakeSpikeOOPSIData(nDataSet);
nYesData             = nDataSet.unit_yes_trial;
nNoData              = nDataSet.unit_no_trial;
stdYesTrial          = std(nYesData, 1)/sqrt(size(nYesData, 1));
stdNoTrial           = std(nNoData, 1)/sqrt(size(nNoData, 1));
nYesData             = mean(nYesData, 1);
nNoData              = mean(nNoData, 1);
maxDff               = max([nYesData, nNoData]);
minDff               = min([nYesData, nNoData]);
diffDff              = maxDff - minDff;        
nYesData             = (nYesData - minDff)/diffDff;
nNoData              = (nNoData - minDff)/diffDff;
stdYesTrial          = stdYesTrial/diffDff;
stdNoTrial           = stdNoTrial/diffDff;
hold on
shadedErrorBar(params.timeSeries, nYesData, stdYesTrial, ...
    {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
shadedErrorBar(params.timeSeries, nNoData, stdNoTrial, ...
    {'k-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
ylabel('Normalized dff (Hz)');
xlabel('Time (s)');
ylim([0 1])
xlim([params.timeSeries(1) params.timeSeries(end)]);
set(gca, 'TickDir', 'out') 
title('C2S model')

setPrint(8*2, 6*2, ['ModelExampleForward_Backward_Trans_' num2str(nNeuron,'%04d') 'nFactor_' num2str(round(nFactor))]) 
setPrint(8*2, 6*2, ['ModelExampleForward_Backward_Trans_' num2str(nNeuron,'%04d') 'nFactor_' num2str(round(nFactor))], 'png')    

close all;

