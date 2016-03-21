%
% Collecting parameters from TCW's file
% Includes:
% 1. Fm
% 2. K
% 3. n
% 4. tau_rise
% 5. tau_decay
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Analysis_Code_v1_TWC 

    addpath('../Func/')
    setDir;
    
    load([TempDir 'TotCell.mat'], 'totCell');
    
    paras       = repmat(struct('cellName',1, 'nRep', 1, 'expression', 'virus',...
                                'CaIndicator', 'GCaMP6f', ...
                                'Fm',1, 'K', 1, 'n', 1,...
                                'tau_r', 1, 'tau_d', 1),length(totCell), 1); 

    numCell     = 0;
    cellToAnalysis = false(length(totCell), 1);
    
    if ~exist([PlotDir 'TWCModelFits/'],'dir')
        mkdir([PlotDir 'TWCModelFits/'])
    end
    

    for nCell   = 1:length(totCell)
        if length(totCell(nCell).CaTime)/(length(totCell(nCell).ephysTime)/1000) == 6
            numCell                      = numCell + 1;
            cellToAnalysis(nCell)        = true;
            paras(numCell).cellName      = totCell(nCell).cellName;
            paras(numCell).nRep          = totCell(nCell).nRep;
            paras(numCell).expression    = totCell(nCell).expression;
            paras(numCell).CaIndicator   = totCell(nCell).CaIndicator;
            para.t_frame                 = totCell(nCell).CaTime;
            para.t_ephys                 = totCell(nCell).ephysTime;
            para.fneuropil               = totCell(nCell).fNeuropil;
            para.fmean                   = totCell(nCell).fROI;
            para.filt                    = totCell(nCell).filteredEphys;
            para.peak                    = totCell(nCell).detectedSpikes;
%             dff                          = get_baseline_corr_dff(para);
            dff                          = get_baseline(para);
            spk                          = {para.t_ephys(para.peak)};
            para_start                   = [20   20.4788    1.1856    1    0.2107];
            para_final                   = gcamp6_model_4para_new(spk, dff, para.t_frame, para_start);
            CaTraces                     = func_getCaTraces_general_new(spk,para.t_frame,para_final);
            paras(numCell).Fm            = para_final(1);
            paras(numCell).K             = para_final(2);
            paras(numCell).n             = para_final(3);
            paras(numCell).tau_d         = para_final(4);
            paras(numCell).tau_r         = para_final(5);
            figure;
            hold on;
            plot(para.t_frame, dff, '-b');
            plot(para.t_frame, CaTraces, '-r')
            plot(spk{1}, ones(size(spk{1}))*7, '+k')
            hold off
            box off
            xlim([para.t_frame(1) para.t_frame(end)])
            ylim([-0.2 8])
            xlabel('Time (s)')
            ylabel('Baseline-corrected F')
            
            setPrint(8,6,[PlotDir 'TWCModelFits/' totCell(nCell).expression '_' ...
                totCell(nCell).cellName '_' num2str(totCell(nCell).nRep,'%02d')],...
                'pdf')
            close all;
        end
    end

    paras       = paras(1:numCell);
    
    save([TempDir 'Paras_TWC_model.mat'], 'paras','cellToAnalysis','-v7.3');    