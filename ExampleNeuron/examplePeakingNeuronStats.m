% shift of peaks

addpath('../Func');
setDir;

load ([TempDatDir 'DataListShuffle.mat']);
nData = 1; % plot raster and psth
load([TempDatDir DataSetList(nData).name '_old.mat'])
spikeDataSet = nDataSet;   
params = DataSetList(nData).params;
%     dynamicalNeuronIndex = [1047 982 972 810 751 742 728 555 401 395 307 264 221 161 137 1041];
%     dynamicalNeuronIndex = dynamicalNeuronIndex([6, 5, 4, 8, 16 ]);
params.Fm = 1;
params.K  = 1;
params.n  = 1;
params.tau_r = 0.060;
params.tau_d = 1.1;
params.intNoise = 0;
params.extNoise = 0;
ext_noise = params.extNoise;
g = @(p,x) p(1) + p(2)./ (1 + 10.^((p(3)-x)*p(4)));
param = [0, 1, 0.2, 3];

numTrial    = 100;

nDelay      = nan(length(spikeDataSet), 2);

for nCell = 1:length(spikeDataSet)
    unit_yes_trial = mean(spikeDataSet(nCell).unit_yes_trial);
    unit_no_trial  = mean(spikeDataSet(nCell).unit_no_trial);

    [max_yes, max_yes_id] = max(unit_yes_trial);
    [max_no,  max_no_id]  = max(unit_no_trial);
    if max_yes > max_no
        max_id  = max_yes_id;
        yesPeak = true;
    else
        max_id  = max_no_id;
        yesPeak = false;
    end

    if max_id <=50 && max_id >=47
        params.tau_d  = 1.7;
        nCellDataSet  = getFakeCaImagingData(spikeDataSet(nCell), params);
        if yesPeak
            unit_trial = mean(nCellDataSet.unit_yes_trial_linear);
            [~, new_max_id] = max(unit_trial);
            nDelay(nCell, 1)   = (new_max_id - max_id)*params.binsize;
        else
            unit_trial = mean(nCellDataSet.unit_no_trial_linear);
            [~, new_max_id] = max(unit_trial);
            nDelay(nCell, 1)   = (new_max_id - max_id)*params.binsize;
        end

        params.tau_d  = 1.1;
        nCellDataSet  = getFakeCaImagingData(spikeDataSet(nCell), params);
        if yesPeak
            unit_trial = mean(nCellDataSet.unit_yes_trial_linear);
            [~, new_max_id] = max(unit_trial);
            nDelay(nCell, 2)   = (new_max_id - max_id)*params.binsize;
        else
            unit_trial = mean(nCellDataSet.unit_no_trial_linear);
            [~, new_max_id] = max(unit_trial);
            nDelay(nCell, 2)   = (new_max_id - max_id)*params.binsize;
        end


    end
end

figure;
plot(nDelay(:,2), nDelay(:,1), 'ok')
hold on
plot([-2 2], [-2 2], '--k')
xlim([-2 2])
ylim([-2 2])
xlabel('S2C')
ylabel('S2C')
box off
set(gca, 'TickDir', 'out')
setPrint(8, 6, [PlotDir 'SingleUnitsImagescWithSort/PeakDelayTaus'])


sigma                         = 0.10 / params.binsize; % 200 ms
filterLength                  = 11;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse); 
modulation                    = nan(length(spikeDataSet), 1);

for nCell = 1:length(spikeDataSet)
    unit_yes_trial = mean(spikeDataSet(nCell).unit_yes_trial);
    unit_no_trial  = mean(spikeDataSet(nCell).unit_no_trial);

    [max_yes, max_yes_id] = max(unit_yes_trial);
    [max_no,  max_no_id]  = max(unit_no_trial);
    if max_yes > max_no
        max_id  = max_yes_id;
        yesPeak = true;
    else
        max_id  = max_no_id;
        yesPeak = false;
    end

    if max_id <=50 && max_id >=47
        if yesPeak
            nUnitData        = spikeDataSet(nCell).unit_yes_trial;
            nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
            y                = mean(nUnitData, 1);
            modulation(nCell)= max_yes/y(max_id);
        else
            nUnitData        = spikeDataSet(nCell).unit_no_trial;
            nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
            y                = mean(nUnitData, 1);
            modulation(nCell)= max_no/y(max_id);
        end
    end
end

figure;
plot(modulation, nDelay(:,1), 'ok')
ylim([-2 2])
xlabel('S2C')
ylabel('S2C')
box off
set(gca, 'TickDir', 'out')
setPrint(8, 6, [PlotDir 'SingleUnitsImagescWithSort/PeakDelayModulation'])