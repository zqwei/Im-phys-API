%
% Comparison between data and TWC model
% Noise in time
% Noise of single AP
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

addpath('../Func/')
setDir;
load([TempDir 'TotCell.mat'], 'totCell');
load([TempDir 'Paras_TWC_model.mat'], 'paras','cellToAnalysis');        
cellToAnalysis    = find(cellToAnalysis);
parasTWCModel     = paras;
preAPLength       = 60;
postAPLength      = 60;
singleAPPattern   = [zeros(preAPLength, 1); 1; zeros(postAPLength, 1)];
singleAPLength    = length(singleAPPattern);
singleAPTime      = preAPLength + 1;
caRate            = 1/60;

dffAll            = {[], []}; % df/f data
modelDffAll       = {[], []}; % df/f from model
modelCaTraceAll   = {[], []}; % linear ca++ signal 
                              % modeled df/f is a nonlinear function of 
                              % linear ca++ signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot of single AP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    subplot(1, length(dffAll), nData);
    plot(caRate * ((1:singleAPLength) - singleAPTime), dff, '-', 'linewid', 0.5, 'color', [0.5 0.5 0.5]);
    hold on;
    plot(caRate * ((1:singleAPLength) - singleAPTime), mean(dff), '-k', 'linewid', 1.5);
    hold off;
    xlim([caRate*-(singleAPTime-1) caRate*(singleAPLength-singleAPTime)]);
    ylabel('df/f')
    xlabel('Time (s)')
    title(caTitle{nData})
end

% for nData            = 1:length(dffAll)
%     dff              = dffAll{nData};
%     modelDff         = modelDffAll{nData};
%     noiseTempCorr    = bsxfun(@times, dff - modelDff, dff(:, singleAPTime) - modelDff(:, singleAPTime));
%     subplot(1, length(dffAll), nData);
%     plot(caRate * ((1:singleAPLength) - singleAPTime), noiseTempCorr(:, :), '-', 'linewid', 0.5, 'color', [0.5 0.5 0.5]);
%     hold on;
%     plot(caRate * ((1:singleAPLength) - singleAPTime), mean(noiseTempCorr), '-k', 'linewid', 1.5);
%     hold off;
%     xlim([caRate*-(singleAPTime-1) caRate*(singleAPLength-singleAPTime)]);
%     ylabel('Temp. Cov.')
%     xlabel('Time (s)')
%     title(caTitle{nData})
%     
% %     xbins       = 0:0.2:12;
% %     [~, nGroup] = histc(modelCaTraceAll{1}(:), xbins);
% %     groupErrs   = grpstats(noiseTempCorr(:), nGroup);
% %     plot(xbins(unique(nGroup)+1), groupErrs,'ob')
% %     xlim([0 12])
% %     box off
% %     xlabel('Modeled [Ca^{++}]')
% %     ylabel('Temp. Cov.')
% %     title(caTitle{nData})
% end
% 
% setPrint(16, 6, [PlotDir 'TWCModelEstimationErrorTemporalCorrelationOneAP'],'pdf')