%
% Example neurons for neuronal trace from S2C model and C2S model
%
%


addpath('../Func');
setDirV1Cells;
load([TempDatDir 'DatasetListCells.mat'], 'totCell');

if ~exist([PlotDir 'ModelCellFits'],'dir')
    mkdir([PlotDir 'ModelCellFits'])
end

paras = repmat(struct('cellName',1, 'nRep', 1, 'expression', 'virus',...
                        'CaIndicator', 'GCaMP6f', 'n_fast', 1, 'n_smc', 1, ...
                        'Fm',1, 'K', 1, 'n', 1, 'tau_r', 1, 'tau_d', 1),length(totCell), 1); 
                    
                    
                    
for nCell   = 1:length(totCell)  
    %% S2C model
    
    spk          = totCell(nCell).spk;
    dff          = totCell(nCell).dff;
    para_start   = [20   20.4788    1.1856    1    0.2107];
    t_frame      = totCell(nCell).CaTime;
    
    %%%%%%%% linear model
    % para_start                   = [20   0    1    0.2107];
    % para_final                   = gcamp6_linear_model({spk}, dff, para.t_frame, para_start);
    % CaTraces                     = func_getCaTraces_linear({spk},para.t_frame,para_final);
    
    %%%%%%%% quadratic model
    % para_start                   = [20  0  0  1  0.2107];
    % para_final                   = gcamp6_quadratic_model({spk}, dff, para.t_frame, para_start);
    % CaTraces                     = func_getCaTraces_quadratic({spk},para.t_frame,para_final);
    
    % Hill model
    para_final   = gcamp6_model_4para_new({spk}, dff, t_frame, para_start);
    fitCaTraces  = func_getCaTraces_general_new({spk}, t_frame,para_final);
    paras(nCell).cellName      = totCell(nCell).cellName;
    paras(nCell).nRep          = totCell(nCell).nRep;
    paras(nCell).expression    = totCell(nCell).expression;
    paras(nCell).CaIndicator   = totCell(nCell).CaIndicator;
    paras(nCell).Fm            = para_final(1);
    paras(nCell).K             = para_final(2);
    paras(nCell).n             = para_final(3);
    paras(nCell).tau_d         = para_final(4);
    paras(nCell).tau_r         = para_final(5);
    
    %% C2S model
    
    V.fast_iter_max    = min(1000, length(t_frame));
    V.smc_iter_max     = 1000;
    V.dt               = t_frame(end)/length(t_frame);
    V.preprocess       = 1; % high-pass filter (increase fitting speed)
    V.T                = length(t_frame);
    [fast, smc]        = runOOPSI(dff', V);
    n_fast             = fast.n/max(fast.n);
    paras(nCell).n_fast = n_fast;
    if ~isempty(smc)
        paras(nCell).n_smc  = smc.E.nbar;
    else
        paras(nCell).n_smc  = n_fast;
    end

            
            
    
    %% plots
    normalized_dff     = (dff - min(dff))/(max(dff)-min(dff));
    normalized_fitCaTraces = (fitCaTraces - min(dff))/(max(dff)-min(dff));
    n_spk              = hist(spk, t_frame);
    
    figure;
    
    subplot(1, 3, 1)    
    hold on;
    plot(para.t_frame, normalized_dff, '-k', 'linewid', 1);
    plot(spk{1}, ones(size(spk{1}))*1.05, 'o', 'color', [0.3 0.3 0.3], 'linewid', 1)
    hold off
    axis off
    xlim([para.t_frame(1) para.t_frame(end)])
    ylim([0 1.1])
    ylabel('original data')

    subplot(1, 3, 2)    
    hold on;
    plot(para.t_frame, normalized_dff, '-k', 'linewid', 1);
    plot(para.t_frame, normalized_fitCaTraces, '-', 'color', [0.3 0.3 0.3], 'linewid', 1)
    hold off
    axis off
    xlim([para.t_frame(1) para.t_frame(end)])
    ylabel('S2C model')
    
    subplot(1, 3, 3)    
    hold on;
    stem(para.t_frame, n_spk, '-k', 'linewid', 1);
    stem(para.t_frame, n_fast, '-', 'linewid', 1)
    stem(para.t_frame, n_smc, '-', 'linewid', 1)
    hold off
    axis off
    xlim([para.t_frame(1) para.t_frame(end)])
    ylabel('C2S model')
    
    setPrint(8,6,[PlotDir 'ModelCellFits/' totCell(nCell).expression '_' ...
        totCell(nCell).cellName '_' num2str(totCell(nCell).nRep,'%02d')])
end

save([TempDatDir 'ParamsFitCells.mat'], 'paras');

close all;