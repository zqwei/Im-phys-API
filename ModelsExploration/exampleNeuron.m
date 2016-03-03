%
% Example neurons for neuronal trace from S2C model and C2S model
%
%


addpath('../Func');
setDirV1Cells;
load([TempDatDir 'DataListCells.mat'], 'totCell');

if ~exist([PlotDir 'ModelExampleCellFits'],'dir')
    mkdir([PlotDir 'ModelExampleCellFits'])
end

% paras = repmat(struct('cellName',1, 'nRep', 1, 'expression', 'virus',...
%                         'CaIndicator', 'GCaMP6f', 'n_fast', 1, 'n_smc', 1, ...
%                         'Fm',1, 'K', 1, 'n', 1, 'tau_r', 1, 'tau_d', 1),length(totCell), 1); 
%                     
%                     
%                     
% for nCell   = 1:length(totCell)  
%     %% S2C model
%     tic
%     
%     spk          = totCell(nCell).spk;
%     dff          = totCell(nCell).dff;
%     if isa(dff, 'single'); dff = double(dff); end
%     para_start   = [20   20.4788    1.1856    1    0.2107];
%     t_frame      = totCell(nCell).CaTime;
%     
%     %%%%%%%% linear model
%     % para_start                   = [20   0    1    0.2107];
%     % para_final                   = gcamp6_linear_model({spk}, dff, para.t_frame, para_start);
%     % CaTraces                     = func_getCaTraces_linear({spk},para.t_frame,para_final);
%     
%     %%%%%%%% quadratic model
%     % para_start                   = [20  0  0  1  0.2107];
%     % para_final                   = gcamp6_quadratic_model({spk}, dff, para.t_frame, para_start);
%     % CaTraces                     = func_getCaTraces_quadratic({spk},para.t_frame,para_final);
%     
%     % Hill model
%     para_final   = gcamp6_model_4para_new({spk}, dff, t_frame, para_start);
%     fitCaTraces  = func_getCaTraces_general_new({spk}, t_frame,para_final);
%     paras(nCell).cellName      = totCell(nCell).cellName;
%     paras(nCell).nRep          = totCell(nCell).nRep;
%     paras(nCell).expression    = totCell(nCell).expression;
%     paras(nCell).CaIndicator   = totCell(nCell).CaIndicator;
%     paras(nCell).Fm            = para_final(1);
%     paras(nCell).K             = para_final(2);
%     paras(nCell).n             = para_final(3);
%     paras(nCell).tau_d         = para_final(4);
%     paras(nCell).tau_r         = para_final(5);
%  
%     %% C2S model
%     
%     V.fast_iter_max    = 1000;
%     V.smc_iter_max     = 1000;
%     V.dt               = t_frame(end)/length(t_frame);
%     V.preprocess       = 1; % high-pass filter (increase fitting speed)
%     V.T                = length(t_frame);
%     [fast, smc]        = runOOPSI(dff', V);
%     n_fast             = fast.n/max(fast.n);
%     paras(nCell).n_fast = n_fast;
%     if ~isempty(smc)
%         n_smc  = smc.E.nbar;
%     else
%         n_smc  = n_fast;
%     end
%     n_smc  = n_smc/max(n_smc);
%     paras(nCell).n_smc  = n_smc;
% 
%             
% end
% 
% save([TempDatDir 'ParamsFitCells.mat'], 'paras');



%% This code is only for plots, where computation is done by server

load([TempDatDir 'ParamsFitCells.mat'], 'paras');

for nCell   = 1:length(totCell)  
    
    spk          = totCell(nCell).spk;
    dff          = totCell(nCell).dff;
    if isa(dff, 'single'); dff = double(dff); end
    t_frame      = totCell(nCell).CaTime;    
    para_final   = [paras(nCell).Fm paras(nCell).K paras(nCell).n paras(nCell).tau_d paras(nCell).tau_r];
    fitCaTraces  = func_getCaTraces_general_new({spk}, t_frame,para_final);

    normalized_dff     = (dff - min(dff))/(max(dff)-min(dff));
    normalized_fitCaTraces = (fitCaTraces - min(dff))/(max(dff)-min(dff));
    n_spk              = hist(spk, t_frame);
    n_spk              = n_spk/max(n_spk);
    
    n_fast             = paras(nCell).n_fast;
    n_smc              = paras(nCell).n_smc; 
    
    figure;
    
    subplot(3, 1, 1)    
    hold on;
    plot(t_frame, normalized_dff, '-k', 'linewid', 1);
    plot(spk, ones(size(spk))*1.05, 'k+', 'linewid', 1)
    hold off
    axis off
    % xlim([t_frame(1) t_frame(end)])
    xlim([0 30])
    ylim([0 1.4])
    title('original data')

    subplot(3, 1, 2)    
    hold on;
    plot(t_frame, normalized_dff, '-k', 'linewid', 1);
    plot(t_frame, normalized_fitCaTraces, '-r', 'linewid', 1)
    plot(3, 1.1, 's', 'color',  'k', 'markerfacecolor', 'k')
    text(3.3, 1.1, 'original DF/F', 'color', 'k')
    plot(13, 1.1, 's', 'color',  'r', 'markerfacecolor', 'r')
    text(13.3, 1.1, 'fit', 'color', 'k')
    hold off
    axis off
    % xlim([t_frame(1) t_frame(end)])
    xlim([0 30])
    ylim([0 1.4])
    title('S2C model')
    
    subplot(3, 1, 3)    
    hold on;
    area(t_frame, n_spk, 'edgecolor', 'k', 'facecolor','k'); alpha(0.5);
    area(t_frame, n_fast, 'edgecolor', [0.4660    0.6740    0.1880], 'facecolor', [0.4660    0.6740    0.1880]); alpha(0.5);
    area(t_frame, n_smc, 'edgecolor', [0.6350    0.0780    0.1840], 'facecolor', [0.6350    0.0780    0.1840]); alpha(0.5);
    plot(3, 1.1, 's', 'color',  'k', 'markerfacecolor', 'k')
    text(3.3, 1.1, 'original firing spikes', 'color', 'k')
    plot(13, 1.1, 's', 'color',  [0.4660    0.6740    0.1880], 'markerfacecolor', [0.4660    0.6740    0.1880])
    text(13.3, 1.1, 'fast fit', 'color', 'k')
    plot(23, 1.1, 's', 'color',  [0.6350    0.0780    0.1840], 'markerfacecolor', [0.6350    0.0780    0.1840])
    text(23.3, 1.1, 'smc fit', 'color', 'k')
    ylim([0 1.4])
    hold off
    axis off
    % xlim([t_frame(1) t_frame(end)])
    xlim([0 30])
    title('C2S model')
    
    setPrint(18, 10,[PlotDir 'ModelExampleCellFits/' totCell(nCell).expression '_' ...
        totCell(nCell).cellName '_' num2str(totCell(nCell).nRep,'%02d')])
end

close all;