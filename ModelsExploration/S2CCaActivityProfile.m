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

function S2CCaActivityProfile
    addpath('../Func/')
    setDir;
    load ([TempDatDir 'DataListShuffle.mat']);
    nData = 1;
    load([TempDatDir DataSetList(nData).name '.mat'])
    if ~exist([PlotDir 'ModelExampleSpikesWithDiffParams'],'dir')
        mkdir([PlotDir 'ModelExampleSpikesWithDiffParams'])
    end
    Fm              = 21.3240;
    K               = 13.9248;
    n               = 1.7531;
    tau_rise        = 0.0728;
    tau_decay       = 1.4551;
    intNoise        = 1.5;
    extNoise        = 0;
    for nNeuron                = 1:length(nDataSet)
        spikeDataSet           = nDataSet(nNeuron);
        params                 = DataSetList(nData).params;
        params.Fm              = Fm;
        params.K               = K;
        params.n               = n;
        params.tau_r           = tau_rise;
        params.tau_d           = tau_decay;
        params.intNoise        = intNoise;
        params.extNoise        = extNoise;   
        figure;   
        spkRaster(spikeDataSet, params);
        spkPSTH(spikeDataSet, params);    
        nKeys   = {'tau_r', 'tau_d', 'intNoise', 'n', 'K'};
        nTitles = {'S2C-linear \tau_{r}', 'S2C-linear \tau_{d}', 'S2C-noise \sigma_{int}', 'S2C-nonlinear n', 'S2C-nonlinear K'};
        nPlots  = [5 8 3 6 9];    
        for nParams = 1:length(nKeys)
            plotChangeParams(spikeDataSet, params, nKeys{nParams}, nTitles{nParams}, nPlots(nParams));
        end
        plotChangeBaseline(spikeDataSet, params);   
        setPrint(8*3, 6*3, [PlotDir 'ModelExampleSpikesWithDiffParams/S2CModelNeuronIndex_' num2str(nNeuron,'%04d')])    
        close all;
    end
end

function spkRaster(spikeDataSet, params)
    spkTimes{1}    = spikeDataSet.unit_yes_trial_spk_time;
    spkTimes{2}    = spikeDataSet.unit_no_trial_spk_time;
    oldBlue         = [     0         0       0.7];
    oldRed          = [    0.7        0         0];
    color_index    = [oldRed; oldBlue];    
    subplot(3, 3, 1)
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
end

function spkPSTH(spikeDataSet, params)
    oldBlue         = [     0         0       0.7];
    oldRed          = [    0.7        0         0];
    color_index     = [oldRed; oldBlue];
    subplot(3, 3, 4)
    hold on;
    nUnitData       = spikeDataSet.unit_yes_trial;
    yesTrialRate    = mean(nUnitData, 1);
    stdYesTrial     = std(nUnitData, 1)/sqrt(size(nUnitData, 1));
    nUnitData       = spikeDataSet.unit_no_trial;
    noTrialRate     = mean(nUnitData, 1);
    stdNoTrial      = std(nUnitData, 1)/sqrt(size(nUnitData, 1));
    maxRate         = max([yesTrialRate, noTrialRate]);  
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
end

function plotChangeParams(spikeDataSet, params, nKey, nTitle, nPlot)
    nValue          = params.(nKey);
    
    oldBlue         = [     0         0       0.7];
    newBlue         = [     0    0.4470    0.7410];
    oldRed          = [    0.7        0         0];
    newRed          = [0.6350    0.0780    0.1840];
    subplot(3, 3, nPlot);
    hold on;
    for nFactor              = 0.5:0.5:1.5
        params.(nKey)        = nValue * nFactor;
        nDataSet             = getFakeCaImagingData(spikeDataSet, params);
        nYesData             = nDataSet.unit_yes_trial;
        nNoData              = nDataSet.unit_no_trial;
        nYesData             = mean(nYesData, 1);
        nNoData              = mean(nNoData, 1);
        maxDff               = max([nYesData, nNoData]);
        minDff               = min([nYesData, nNoData]);
        diffDff              = maxDff - minDff;        
        nYesData             = (nYesData - minDff)/diffDff;
        nNoData              = (nNoData - minDff)/diffDff;
        h = plot(params.timeSeries, nYesData, 'linewid', 1.0);
        switch nFactor
            case 0.5
                h.Color = newRed;
                h.LineStyle = '--';
            case 1.0
                h.Color = oldRed;
                h.LineStyle = '-';
            case 1.5
                h.Color = newRed;
                h.LineStyle = '-';
        end
        h = plot(params.timeSeries, nNoData, 'linewid', 1.0);
        switch nFactor
            case 0.5
                h.Color = newBlue;
                h.LineStyle = '--';
            case 1.0
                h.Color = oldBlue;
                h.LineStyle = '-';
            case 1.5
                h.Color = newBlue;
                h.LineStyle = '-';
        end
    end
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    hold off;
    ylabel('normalized DF/F');
    ylim([0 1])
    xlabel('Time (s)');
    title(nTitle)
end

% function plotChangeBaseline_v1(spikeDataSet, params) %#ok<DEFNU>
%     nDataSet = [spikeDataSet; spikeDataSet; spikeDataSet];
%     nKeys   = {'tau_r', 'tau_d', 'Fm', 'n', 'K'};
%     for n   = 1:length(nKeys)
%         nKey          = nKeys{n};
%         params.(nKey) = params.(nKey) * ones(3, 1);
%     end
%     spkTimes{1}    = spikeDataSet.unit_yes_trial_spk_time;
%     spkTimes{2}    = spikeDataSet.unit_no_trial_spk_time;
%     newSpkTimes1   = spkTimes;
%     newSpkTimes2   = spkTimes;
%     for nType        = 1:length(spkTimes)
%         for nTrial   = 1:length(spkTimes{nType})
%             spk      = spkTimes{nType}{nTrial};
%             if length(spk)>1
%                 DiffSpk  = diff(spk);
%                 addSpk1  = DiffSpk.*rand(size(DiffSpk)) + spk(1:end-1);
%                 newSpkTimes1{nType}{nTrial} = [spk; addSpk1];
%                 addSpk2  = DiffSpk.*rand(size(DiffSpk)) + spk(1:end-1);
%                 newSpkTimes2{nType}{nTrial} = [spk; addSpk1; addSpk2];
%             end
%         end
%     end
%     subplot(3, 3, 2);
%     hold on;
%     oldBlue         = [     0         0       0.7];
%     newBlue         = [     0    0.4470    0.7410];
%     oldRed          = [    0.7        0         0];
%     newRed          = [0.6350    0.0780    0.1840];
%     nDataSet(2).unit_yes_trial_spk_time = newSpkTimes1{1};
%     nDataSet(3).unit_yes_trial_spk_time = newSpkTimes2{1};
%     nDataSet(2).unit_no_trial_spk_time  = newSpkTimes1{2};
%     nDataSet(3).unit_no_trial_spk_time  = newSpkTimes2{2};
%     mDataSet         = getFakeCaImagingData(nDataSet, params);    
%     for mData        = 1:length(mDataSet)
%         nYesData             = mDataSet(mData).unit_yes_trial;
%         nNoData              = mDataSet(mData).unit_no_trial;
%         nYesData             = mean(nYesData, 1);
%         nNoData              = mean(nNoData, 1);        
%         maxDff               = max([nYesData, nNoData]);
%         minDff               = min([nYesData, nNoData]);
%         diffDff              = maxDff - minDff;        
%         nYesData             = (nYesData - minDff)/diffDff;
%         nNoData              = (nNoData - minDff)/diffDff;        
%         h = plot(params.timeSeries, nYesData, 'linewid', 1.0);
%         switch mData
%             case 2
%                 h.Color = newRed;
%                 h.LineStyle = '--';
%             case 1
%                 h.Color = oldRed;
%                 h.LineStyle = '-';
%             case 3
%                 h.Color = newRed;
%                 h.LineStyle = '-';
%         end
%         h = plot(params.timeSeries, nNoData, 'linewid', 1.0);
%         switch mData
%             case 2
%                 h.Color = newBlue;
%                 h.LineStyle = '--';
%             case 1
%                 h.Color = oldBlue;
%                 h.LineStyle = '-';
%             case 3
%                 h.Color = newBlue;
%                 h.LineStyle = '-';
%         end
%     end
%     gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
%     xlim([params.timeSeries(1) params.timeSeries(end)]);
%     hold off;
%     ylabel('normalized DF/F');
%     ylim([0 1])
%     xlabel('Time (s)');
%     title('baseline')
% end

function plotChangeBaseline(spikeDataSet, params)   
    timeBins = params.timeSeries;
    binSize  = params.binsize;
    nDataSet = [spikeDataSet; spikeDataSet; spikeDataSet];
    nKeys    = {'tau_r', 'tau_d', 'Fm', 'n', 'K'};
    for n    = 1:length(nKeys)
        nKey          = nKeys{n};
        params.(nKey) = params.(nKey) * ones(3, 1);
    end

    nUnitData       = spikeDataSet.unit_yes_trial;
    yesTrialRate    = mean(nUnitData, 1);
    numYesTrial     = size(nUnitData, 1);
    nUnitData       = spikeDataSet.unit_no_trial;
    noTrialRate     = mean(nUnitData, 1);
    numNoTrial      = size(nUnitData, 1);

    for nFactor     = 2:3
        [spkTime, spkCounts] = NHpoisson(yesTrialRate * nFactor, timeBins, binSize, numYesTrial);
        nDataSet(nFactor).unit_yes_trial          = spkCounts;
        nDataSet(nFactor).unit_yes_trial_spk_time = spkTime;
        [spkTime, spkCounts] = NHpoisson(noTrialRate * nFactor, timeBins, binSize, numNoTrial);
        nDataSet(nFactor).unit_no_trial           = spkCounts;
        nDataSet(nFactor).unit_no_trial_spk_time  = spkTime;
    end

    subplot(3, 3, 7)
    yesColors      = [    0.7        0         0
                    0.6350    0.0780    0.1840
                    0.6350    0.0780    0.1840];
    noColors       = [     0         0       0.7
                    0    0.4470    0.7410
                    0    0.4470    0.7410];
    lineStyles    = {'-', '--', '-'};
    hold on;
    for nFactor     = 1:3
        nUnitData       = nDataSet(nFactor).unit_yes_trial;
        yesTrialRate    = mean(nUnitData, 1);
        stdYesTrial     = std(nUnitData, 1)/sqrt(size(nUnitData, 1));
        nUnitData       = nDataSet(nFactor).unit_no_trial;
        noTrialRate     = mean(nUnitData, 1);
        stdNoTrial      = std(nUnitData, 1)/sqrt(size(nUnitData, 1));
        maxRate         = max([yesTrialRate, noTrialRate]);  
        shadedErrorBar(params.timeSeries, yesTrialRate/maxRate, stdYesTrial/maxRate, ...
            {lineStyles{nFactor}, 'linewid', 1.0, 'color', yesColors(nFactor, :)}, 0.5);
        shadedErrorBar(params.timeSeries, noTrialRate/maxRate, stdNoTrial/maxRate, ...
            {lineStyles{nFactor}, 'linewid', 1.0, 'color', noColors(nFactor, :)}, 0.5);
    end
    hold off;
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    hold off;
    ylabel('normalized rate');
    ylim([0 1])
    xlabel('Time (s)');
    title('baseline')


    subplot(3, 3, 2);
    hold on;
    oldBlue         = [     0         0       0.7];
    newBlue         = [     0    0.4470    0.7410];
    oldRed          = [    0.7        0         0];
    newRed          = [0.6350    0.0780    0.1840];
    mDataSet                 = getFakeCaImagingData(nDataSet, params);    
    for mData                = 1:length(mDataSet)
        nYesData             = mDataSet(mData).unit_yes_trial;
        nNoData              = mDataSet(mData).unit_no_trial;
        nYesData             = mean(nYesData, 1);
        nNoData              = mean(nNoData, 1);        
        maxDff               = max([nYesData, nNoData]);
        minDff               = min([nYesData, nNoData]);
        diffDff              = maxDff - minDff;        
        nYesData             = (nYesData - minDff)/diffDff;
        nNoData              = (nNoData - minDff)/diffDff;        
        h = plot(params.timeSeries, nYesData, 'linewid', 1.0);
        switch mData
            case 2
                h.Color = newRed;
                h.LineStyle = '--';
            case 1
                h.Color = oldRed;
                h.LineStyle = '-';
            case 3
                h.Color = newRed;
                h.LineStyle = '-';
        end
        h = plot(params.timeSeries, nNoData, 'linewid', 1.0);
        switch mData
            case 2
                h.Color = newBlue;
                h.LineStyle = '--';
            case 1
                h.Color = oldBlue;
                h.LineStyle = '-';
            case 3
                h.Color = newBlue;
                h.LineStyle = '-';
        end
    end
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    hold off;
    ylabel('normalized DF/F');
    ylim([0 1])
    xlabel('Time (s)');
    title('baseline')
end

function [spkTime, spkCounts] = NHpoisson(rate, timeBins, binSize, numTrial)
    % algorithm from http://transp-or.epfl.ch/courses/OptSim2012/slides/05b-poisson.pdf
    dt         = 0.001; % sampling rate = 1000 Hz
    t_start    = timeBins(1);
    t_end      = timeBins(end);
    maxRate    = max(rate);
    spkTime    = cell(numTrial, 1);
    spkCounts  = zeros(numTrial, length(timeBins));
    for nTrial = 1:numTrial
        kSpk     = 0;
        while kSpk == 0
            t        = t_start;
            kSpk     = 0;
            spkTrial = nan(ceil((t_end-t_start)/dt), 1);
            while t<=t_end
                r    = rand();
                t    = t - log(r)/maxRate;
                s    = rand();
                lambda_t = rate(sum(t>timeBins));
                if s<= lambda_t/maxRate
                    kSpk= kSpk+1;
                    spkTrial(kSpk) = t;
                end
            end
        end
        spkTrial             = spkTrial(~isnan(spkTrial));
        spkTime{nTrial}      = spkTrial;
        spkCounts(nTrial, :) = hist(spkTrial,timeBins)/binSize;
    end
end