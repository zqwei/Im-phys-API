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

function S2CActivityProfileExample
    addpath('../Func/')
    setDir;
    load ([TempDatDir 'ParamsFitCells_S2CModel_Sim.mat'], 'params');
    paramsSim = params; %#ok<NODEF>
    load ([TempDatDir 'DataListShuffle.mat']);
    nData = 1;
    load([TempDatDir DataSetList(nData).name '.mat'])
    if ~exist([PlotDir 'ModelCellFits'],'dir')
        mkdir([PlotDir 'ModelCellFits'])
    end
    params     = DataSetList(nData).params;
    numT       = size(nDataSet(1).unit_yes_trial,2);
    actMat     = nan(length(nDataSet),numT*2);
    for nUnit  = 1: length(nDataSet)
        actMat(nUnit, 1:numT)     = mean(nDataSet(nUnit).unit_yes_trial);
        actMat(nUnit, numT+1:end) = mean(nDataSet(nUnit).unit_no_trial);
    end
    actMat     = mean(actMat, 1);    
    actYes     = actMat(2:numT-1);
    actYes     = min(actYes) + (actYes - min(actYes))/(max(actYes)-min(actYes)) * 20;
    actNo      = actMat(numT+2:end-1);
    actNo      = min(actNo) + (actNo - min(actNo))/(max(actNo)-min(actNo)) * 20;
    actMat     = [actYes; actNo];
    actName    = {'ipsi', 'contra'};
    timeBins   = params.timeSeries(2:end-1);
    numTrial   = 100;
    params.Fm  = 21.3240;
    params.extNoise = 0;
    params.intNoise = 0.5;
    binSize    = params.binsize;
    
    groupNames      = {'6f_AAV', '6s_AAV', '6f_Tsg', '6s_Tsg'};
    colorSet        = cbrewer('div', 'Spectral', 10, 'cubic');
    nKeys           = {'tau_r', 'tau_d', 'n', 'K'};
    nKeyNames       = {'Change of \tau_{rise}', 'Change of \tau_{decay}', 'Change of n', 'Change of k', 'Change of firing rate'};
    yLabelNames     = {'\tau_{r} (ms)', '\tau_{d} (ms)', 'n', 'k', 'max rate ratio'};
    nScales         = [1000, 1000, 1, 1];

%     nKeys           = {'tau_d'};
%     nKeyNames       = {'Change of \tau_{decay}', 'Change of firing rate'};
%     yLabelNames     = {'\tau_{d} (ms)', 'max rate ratio'};
%     nScales         = [1000];
    
    for nGroup      = [2 4]        
        params.K               = paramsSim(1, nGroup).K;
        params.n               = paramsSim(1, nGroup).n;
        params.tau_r           = paramsSim(1, nGroup).tau_r;
        params.tau_d           = paramsSim(1, nGroup).tau_d;
        
        for nAct = 1:size(actMat, 1)   
            figure;
            for nKey     = 1:length(nKeys)
                tparams  = params;
                % change of decay time constant
                [spkTime, spkCounts] = NHpoisson(actMat(nAct, :), timeBins, binSize, numTrial);            
                numStep            = 10;
                spkCountMat        = nan(numStep, size(spkCounts, 2));
                nFactor            = linspace(paramsSim(2, nGroup).(nKeys{nKey}), paramsSim(3, nGroup).(nKeys{nKey}), numStep);
                spkTimes           = cell(numStep, numTrial);    
                DFFMat             = nan(numStep, size(spkCounts, 2));
                for nn             = 1:10
                    tparams.(nKeys{nKey})   = nFactor(nn);
                    spkTimes(nn, :)= spkTime;
                    spkCountMat(nn, :) = (mean(spkCounts, 1)-min(mean(spkCounts, 1)))./(max(mean(spkCounts, 1))--min(mean(spkCounts, 1)));
                    DFFMat(nn, :)  = S2CModel(spkTime, spkCounts, tparams);
                end                     
                subplot(3, length(nKeyNames), nKey)
                title(nKeyNames{nKey})
                spkRaster(spkTimes, tparams, colorSet)            
                subplot(3, length(nKeyNames), nKey + length(nKeyNames))
                title('normalized firing')
                spkHistMat(spkCountMat, tparams, yLabelNames{nKey}, nFactor*nScales(nKey))
                subplot(3, length(nKeyNames), nKey + length(nKeyNames)*2)
                title('normalized dff')
                spkHistMat(DFFMat, tparams, yLabelNames{nKey}, nFactor*nScales(nKey))
            end
            % change of spk rate
            numStep            = 10;
            spkCountMat        = nan(numStep, size(spkCounts, 2));
            spkTimes           = cell(numStep, numTrial);    
            DFFMat             = nan(numStep, size(spkCounts, 2));
            for nn             = 1:numStep
                [spkTime, spkCounts] = NHpoisson(actMat(nAct, :)*((nn-1)*0.4+1), timeBins, binSize, numTrial);            
                spkTimes(nn, :)= spkTime;
                spkCountMat(nn, :) = (mean(spkCounts, 1)-min(mean(spkCounts, 1)))./(max(mean(spkCounts, 1))--min(mean(spkCounts, 1)));
                DFFMat(nn, :)  = S2CModel(spkTime, spkCounts, params);
            end
            nn = 1:numStep;
            subplot(3, length(nKeyNames), nKey+1)
            title('Change of max firing rate')
            spkRaster(spkTimes, params, colorSet)            
            subplot(3, length(nKeyNames), nKey+1 + length(nKeyNames))
            title('normalized firing')
            spkHistMat(spkCountMat, params, yLabelNames{nKey+1}, (nn-1)*0.4+1)
            subplot(3, length(nKeyNames), nKey+1 + length(nKeyNames)*2)
            title('normalized dff')
            spkHistMat(DFFMat, params, yLabelNames{nKey+1}, (nn-1)*0.4+1)
            setPrint(8*(nKey + 1), 6*3, [PlotDir 'ModelCellFits/S2CActivityProfile_' groupNames{nGroup} '_' actName{nAct} '_' num2str(length(nKeys)) '_keys'], 'tif');  
        end
    end
    close all
end

function spkRaster(spkTime, params, spkColor)
    hold on;
    spkLoc           = 0;
    for nPlot        = 1:size(spkTime, 1)
        for nTrial       = 1:80/size(spkTime, 1)
            spikeTimes   = spkTime{nPlot, nTrial};
            spikeTrial   = ones(length(spikeTimes), 1) * (nTrial + spkLoc);
            plot(spikeTimes, spikeTrial, '.', 'color', spkColor(nPlot, :));
        end
        spkLoc     = spkLoc + nTrial + 3;
    end
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0);
    hold off;
    ylim([1 spkLoc])
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    axis off
end

function spkHistMat(spkCountMat, params, ylabelName,yTicks)
    hold on;
    imagesc(params.timeSeries(2:end-1), yTicks, spkCountMat, [0 1])
    axis xy
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel(ylabelName);
    xlabel('Time (s)');
    xlim([params.timeSeries(2) params.timeSeries(end-1)]);
    ylim([yTicks(1) yTicks(end)]);
    if ceil(yTicks(1))~= floor(yTicks(end))
        set(gca, 'yTick', [ceil(yTicks(1)) floor(yTicks(end))])
    else
        set(gca, 'yTick', [ceil(yTicks(1)*10)/10 floor(yTicks(end)*10)/10])
    end
    set(gca, 'TickDir', 'out')
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

function DFF = S2CModel(spkTime, spkCounts, params)
                            
    timePoints         = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries(2:end-1));
    timeSeriesData     = params.timeSeries(2:end-1);
    intNoise           = params.intNoise;
    constNoise         = params.extNoise;
    
    paramsSet      = [params.Fm, params.K, params.n, params.tau_d, params.tau_r, intNoise];
    rMean          = mean(mean(spkCounts(:,timePoints(1):timePoints(2))));
    nDFF           = spikeTimeToImaging(spkTime, timeSeriesData, paramsSet, rMean);
    DFF            = nDFF + randn(size(nDFF))*constNoise;
    DFF            = mean(DFF, 1);
    DFF            = (DFF - min(DFF))/(max(DFF) - min(DFF));
end
