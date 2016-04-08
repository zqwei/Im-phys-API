%
% Temporal correlation of the noise using TWC model
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

% function Analysis_Code_v1_0_3 

    addpath('../Func/')
    setDir;
    load([TempDir 'TotCell.mat'], 'totCell');
    load([TempDir 'Paras_TWC_model.mat'], 'paras','cellToAnalysis');        
    cellToAnalysis    = find(cellToAnalysis);
    
    if ~exist([PlotDir 'TWCModelXCorr/'],'dir')
        mkdir([PlotDir 'TWCModelXCorr/'])
    end

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
        [CaTraces, CaTracesOrg]      = func_getCaTraces_general_new(spk, para.t_frame, para_final);
        
        [acor,lag]                   = xcorr(dff - CaTraces);
        
        plot(lag, acor)
        
        setPrint(8,6,[PlotDir 'TWCModelXCorr/' totCell(nCell).expression '_' ...
                totCell(nCell).cellName '_' num2str(totCell(nCell).nRep,'%02d')],...
                'pdf')
    end