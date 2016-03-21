%
% Comparison between data and TWC model
% Noise in time
% Noise of two APs
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

% function Analysis_Code_v1_0_1 

addpath('../Func/')
setDir;
load([TempDir 'TotCell.mat'], 'totCell');
load([TempDir 'Paras_TWC_model.mat'], 'paras','cellToAnalysis');        
cellToAnalysis    = find(cellToAnalysis);
parasTWCModel     = paras;

twoAPPattern      = [zeros(60, 1); 1; zeros(10, 1); 1; zeros(60, 1)]; % one can change this
singleAPLength    = 131;
singleAPTime      = 11;
caRate            = 1/60;

dffAll            = {[], []};
modelDffAll       = {[], []};
modelCaTraceAll   = {[], []};

for tCell         = 1:length(parasTWCModel)
    nCell                        = cellToAnalysis(tCell);
    para.t_frame                 = totCell(nCell).CaTime;
    para.t_ephys                 = totCell(nCell).ephysTime;
    para.fneuropil               = totCell(nCell).fNeuropil;
    para.fmean                   = totCell(nCell).fROI;
    para.filt                    = totCell(nCell).filteredEphys;
    para.peak                    = totCell(nCell).detectedSpikes;
    dff                          = get_baseline(para);
    spk                          = {para.t_ephys(para.peak)};
    para_final                   = [paras(tCell).Fm, paras(tCell).K, ...
                                    paras(tCell).n, paras(tCell).tau_d, ...
                                    paras(tCell).tau_r];
    [CaTraces, CaTracesOrg]      = func_getCaTraces_general_new(spk, para.t_frame, para_final);
    caPeak                       = zeros(size(dff));
    spikeTimes                   = find(para.peak);
    for nSpike                   = 1:length(spikeTimes);
        nSpikeTime               = para.t_ephys(spikeTimes(nSpike));
        nCaTime                  = sum(nSpikeTime >= para.t_frame);
        caPeak(nCaTime)          = caPeak(nCaTime) + 1;
    end
    
    singleAPIndex                = strfind(caPeak', singleAPPattern');

    if strcmp(parasTWCModel(tCell).CaIndicator,'GCaMP6f')
        for nAP                  = 1:length(singleAPIndex)
            dffAll{1}            = [dffAll{1}; dff(nAP:nAP+singleAPLength-1)'];
            modelDffAll{1}       = [modelDffAll{1}; CaTraces(nAP:nAP+singleAPLength-1)'];
            modelCaTraceAll{1}   = [modelCaTraceAll{1}; CaTracesOrg(nAP:nAP+singleAPLength-1)'];
        end
    else
        for nAP                  = 1:length(singleAPIndex)
            dffAll{2}            = [dffAll{2}; dff(nAP:nAP+singleAPLength-1)'];
            modelDffAll{2}       = [modelDffAll{2}; CaTraces(nAP:nAP+singleAPLength-1)'];
            modelCaTraceAll{2}   = [modelCaTraceAll{2}; CaTracesOrg(nAP:nAP+singleAPLength-1)'];
        end
    end
end

caTitle              = {'GCaMP6f', 'GCaMP6s'};
figure;
for nData            = 1:length(dffAll)
    dff              = dffAll{nData};
    modelDff         = modelDffAll{nData};
    noiseTempCorr    = bsxfun(@times, dff - modelDff, dff(:, singleAPTime) - modelDff(:, singleAPTime));
    subplot(1, length(dffAll), nData);
    plot(caRate * ((1:singleAPLength) - 11), noiseTempCorr(1:3, :), '-', 'linewid', 0.5, 'color', [0.5 0.5 0.5]);
    hold on;
    plot(caRate * ((1:singleAPLength) - 11), mean(noiseTempCorr), '-k', 'linewid', 1.5);
    hold off;
    xlim([caRate*-10 caRate*(singleAPLength-11)]);
    ylabel('Temp. Cov.')
    xlabel('Time (s)')
    title(caTitle{nData})
    
%     xbins       = 0:0.2:12;
%     [~, nGroup] = histc(modelCaTraceAll{1}(:), xbins);
%     groupErrs   = grpstats(noiseTempCorr(:), nGroup);
%     plot(xbins(unique(nGroup)+1), groupErrs,'ob')
%     xlim([0 12])
%     box off
%     xlabel('Modeled [Ca^{++}]')
%     ylabel('Temp. Cov.')
%     title(caTitle{nData})
end

setPrint(16, 6, [PlotDir 'TWCModelEstimationErrorTemporalCorrelation'],'pdf')