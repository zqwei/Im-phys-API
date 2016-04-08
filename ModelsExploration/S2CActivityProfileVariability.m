% 
% This code is a summary of results in S2CCaActivityProfile.m (for single
% neurons)
%
%
% version 1.0:
% using doubling of parameters
%
% version 1.1:
% using 10% percentile and 90% percentile of the extreme data for
% parameters
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function S2CActivityProfileVariability
    addpath('../Func/')
    setDir;
    load ([TempDatDir 'DataListShuffle.mat']);
    load ([TempDatDir 'ParamsFitCells_S2CModel_Sim.mat'], 'params');
    paramsSim = params; %#ok<NODEF>
    clear params;
    nData = 1;
    load([TempDatDir DataSetList(nData).name '.mat'])
    if ~exist([PlotDir 'ModelCellFits'],'dir')
        mkdir([PlotDir 'ModelCellFits'])
    end

    Fm              = 21.3240;
    intNoise        = 1.5;
    extNoise        = 0;
    groupNames      = {'6f_AAV', '6s_AAV', '6f_Tsg', '6s_Tsg'};
    nKeys           = {'tau_r', 'tau_d', 'intNoise', 'n', 'K'};
    for nGroup      = [2 4]
        varNParams      = nan(length(nDataSet), length(nKeys)+2, 2);
        for nNeuron                = 1:length(nDataSet)
            disp(nNeuron)
            spikeDataSet           = nDataSet(nNeuron);
            params                 = DataSetList(nData).params;
            params.Fm              = Fm;
            params.K               = paramsSim(1, nGroup).K;
            params.n               = paramsSim(1, nGroup).n;
            params.tau_r           = paramsSim(1, nGroup).tau_r;
            params.tau_d           = paramsSim(1, nGroup).tau_d;
            params.intNoise        = intNoise;
            params.extNoise        = extNoise;   
            nParamsSim             = paramsSim(:, nGroup);
            for ii                 = 1:length(nParamsSim)
                nParamsSim(ii).intNoise = intNoise * ii/2;
            end                
            for nParams = 1:length(nKeys)
                varNParams(nNeuron, nParams+2, :) = varChangeParams(spikeDataSet, params, nParamsSim, nKeys{nParams});
            end
            [varNParams(nNeuron, 1, :), varNParams(nNeuron, 2, :)] = varChangeBaseline(spikeDataSet, params);
        end    
        save([TempDatDir 'PValueModelVariability' groupNames{nGroup} '.mat'], 'varNParams')
        load([TempDatDir 'PValueModelVariability' groupNames{nGroup} '.mat'], 'varNParams')
        nTitles = {'rate', 'baseline', '\tau_{r}', '\tau_{d}', '\sigma_{int}', 'n', 'K'};
        plotCorrMatrix(varNParams, nTitles)
        setPrint(8*2, 8*2, [PlotDir 'ModelCellFits/S2CVariabilityActivityProfile_' groupNames{nGroup}])
        setPrint(8*2, 8*2, [PlotDir 'ModelCellFits/S2CVariabilityActivityProfile_' groupNames{nGroup}],'png')
    end
end

function plotCorrMatrix(varNParams, nTitles)    
    groupColor = [0.7 0 0; 0 0 0.7];
    figure;
    numNeuron = size(varNParams, 1);
    group     = [true(numNeuron, 1); false(numNeuron, 1)];
    xVar      = [squeeze(varNParams(:, :, 1)); squeeze(varNParams(:, :, 2))]; 
    xVar      = -log(xVar);
    [~, ax, ~] = gplotmatrix(xVar, [], group, groupColor, '..', [], 'off', [], nTitles, nTitles);
    %     xlabel(hax, '-log(P-value)')
    %     ylabel(hax, '-log(P-value)')
    xyLimMin   = 0;
    xyLimMax   = 50;
    
    for nAx   = 1:length(nTitles)
        for mAx = 1:length(nTitles)
            ylim(ax(mAx, nAx), [xyLimMin xyLimMax]);
            xlim(ax(mAx, nAx), [xyLimMin xyLimMax]);            
            % non-diagonal plots
            if nAx ~= mAx
                axes(ax(mAx, nAx)); %#ok<LAXES>
                hold on;
                gridxy (3, 3, 'Color','k','Linestyle','--','linewid', 1.0);
                hold off;
            end            
    %             set(ax(mAx, nAx), 'XTick', [0 log(10^10) log(10^20)])
    %             if ~isempty(get(ax(mAx, nAx), 'XTickLabel'))
    %                 set(ax(mAx, nAx), 'XTickLabel', {'1', '10^{-10}', '10^{-20}'})
    %             end
    % 
    %             set(ax(mAx, nAx), 'YTick', [0 log(10^10) log(10^20)])
    %             if ~isempty(get(ax(mAx, nAx), 'YTickLabel'))
    %                 set(ax(mAx, nAx), 'YTickLabel', {'1', '10^{-10}', '10^{-20}'})
    %             end            
        end       
        mAx    = length(nTitles) + 1;
        cla(ax(mAx, nAx))
        axes(ax(mAx, nAx)); %#ok<LAXES>
        hold on;
        sVar   = -log(squeeze(varNParams(:, nAx, :)));
        for nCol = 1:size(sVar, 2)
            [f, xi] = ksdensity(sVar(:, nCol), linspace(xyLimMin, xyLimMax, 40));
            plot(xi, f, 'color', groupColor(nCol, :), 'linewid', 2.0);
        end
        xlim(ax(mAx, nAx), [xyLimMin xyLimMax]);
        gridxy (3, [], 'Color','k','Linestyle','--','linewid', 1.0);
        hold off;       
    %         set(ax(mAx, nAx), 'XTick', [0 log(10^10) log(10^20)])
    %         if ~isempty(get(ax(mAx, nAx), 'XTickLabel'))
    %             set(ax(mAx, nAx), 'XTickLabel', {'1', '10^{-10}', '10^{-20}'})
    %         end
    %         
    %         set(ax(mAx, nAx), 'YTick', [0 log(10^10) log(10^20)])
    %         if ~isempty(get(ax(mAx, nAx), 'YTickLabel'))
    %             set(ax(mAx, nAx), 'YTickLabel', {'1', '10^{-10}', '10^{-20}'})
    %         end
    end
end

function varNParams = varChangeParams(spikeDataSet, params, paramsSim, nKey)
    YesData         = zeros(3, length(params.timeSeries));
    NoData          = zeros(3, length(params.timeSeries));        
    for nFactor              = 1:3
        params.(nKey)        = paramsSim(nFactor).(nKey);
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
        YesData(nFactor, :) = nYesData;
        NoData(nFactor, :)  = nNoData;
    end
    varNParams(1)            = kruskalwallis(YesData', [], 'off');
    varNParams(2)            = kruskalwallis(NoData', [], 'off');    
end

function [varNParams1, varNParams2] = varChangeBaseline(spikeDataSet, params)
    timeBins = params.timeSeries;
    binSize  = params.binsize;        
    nDataSet = [spikeDataSet; spikeDataSet; spikeDataSet];
    tKeys   = {'tau_r', 'tau_d', 'Fm', 'n', 'K'};
    for n   = 1:length(tKeys)
        nKey          = tKeys{n};
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

    YesData          = zeros(3, length(params.timeSeries));
    NoData           = zeros(3, length(params.timeSeries));    
    for mData        = 1:length(nDataSet)
        nYesData             = nDataSet(mData).unit_yes_trial;
        nNoData              = nDataSet(mData).unit_no_trial;
        nYesData             = mean(nYesData, 1);
        nNoData              = mean(nNoData, 1);        
        maxDff               = max([nYesData, nNoData]);
        nYesData             = nYesData/maxDff;
        nNoData              = nNoData/maxDff;
        YesData(mData, :)    = nYesData;
        NoData(mData, :)     = nNoData;
    end   
    varNParams1(1)            = kruskalwallis(YesData', [], 'off');
    varNParams1(2)            = kruskalwallis(NoData', [], 'off');  
    
    mDataSet         = getFakeCaImagingData(nDataSet, params);    
    YesData          = zeros(3, length(params.timeSeries));
    NoData           = zeros(3, length(params.timeSeries));    
    for mData        = 1:length(mDataSet)
        nYesData             = mDataSet(mData).unit_yes_trial;
        nNoData              = mDataSet(mData).unit_no_trial;
        nYesData             = mean(nYesData, 1);
        nNoData              = mean(nNoData, 1);        
        maxDff               = max([nYesData, nNoData]);
        minDff               = min([nYesData, nNoData]);
        diffDff              = maxDff - minDff;
        nYesData             = (nYesData - minDff)/diffDff;
        nNoData              = (nNoData - minDff)/diffDff;
        YesData(mData, :)    = nYesData;
        NoData(mData, :)     = nNoData;
    end    
    varNParams2(1)            = kruskalwallis(YesData', [], 'off');
    varNParams2(2)            = kruskalwallis(NoData', [], 'off');    
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