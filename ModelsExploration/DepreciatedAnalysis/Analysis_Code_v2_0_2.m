%
% Comparison between linear model and TWC model
% noise vs modeled data
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

    dffAll            = {[], []};
    modelDffAll       = {[], []};

    for tCell         = 1:length(paras)
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
        CaTraces                     = func_getCaTraces_general_new(spk, para.t_frame, para_final);
        
        if strcmp(parasTWCModel(tCell).CaIndicator,'GCaMP6f')
            dffAll{1}                = [dffAll{1}; dff];
            modelDffAll{1}           = [modelDffAll{1}; CaTraces];
        else
            dffAll{2}                = [dffAll{2}; dff];
            modelDffAll{2}           = [modelDffAll{2}; CaTraces];
        end
    end
    
    figure;
    subplot(2,2,1)
    hold on
    plot(modelDffAll{1}, dffAll{1} - modelDffAll{1}, '.k')
    plot([0 12], [0 0], '--r')
    xlabel('modeled df/f')
    ylabel('estimation error')
    title('GCaMP6f')
    box off
    xlim([0 12])
    subplot(2,2,2)
    hold on
    plot(modelDffAll{2}, dffAll{2} - modelDffAll{2}, '.k')    
    plot([0 12], [0 0], '--r')
    xlabel('modeled df/f')
    ylabel('estimation error')
    title('GCaMP6s')
    xlim([0 12])
    box off
    
    subplot(2,2,3)
    xbins       = -0.2:0.2:12;
    [~, nGroup] = histc(modelDffAll{1}, xbins);
    groupErrs   = grpstats((dffAll{1} - modelDffAll{1}).^2, nGroup);
    plot(xbins(unique(nGroup)), sqrt(groupErrs),'ob')
    xlim([0 12])
    box off
    xlabel('modeled df/f')
    ylabel('std of estimation error')
    title('GCaMP6f')

    subplot(2,2,4)
    xbins       = -0.2:0.2:12;
    [~, nGroup] = histc(modelDffAll{2}, xbins);
    groupErrs   = grpstats((dffAll{2} - modelDffAll{2}).^2, nGroup);
    plot(xbins(unique(nGroup)), sqrt(groupErrs),'ob')
    xlim([0 12])
    box off
    xlabel('modeled df/f')
    ylabel('std of estimation error')
    title('GCaMP6s')
    
    setPrint(16, 12, [PlotDir 'TWCModelEstimationError'],'pdf')