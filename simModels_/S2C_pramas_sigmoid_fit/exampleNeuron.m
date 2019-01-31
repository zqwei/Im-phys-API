%
% Example neurons for neuronal trace from S2C model and C2S model
% This code is only for plots, where computation is done by server
%
%

addpath('../Func');
setDirV1Cells;
load([TempDatDir 'DataListCells.mat'], 'totCell');

if ~exist([PlotDir 'ModelExampleCellFits'],'dir')
    mkdir([PlotDir 'ModelExampleCellFits'])
end

load([TempDatDir 'ParamsFitCells_S2CModel_sigmoid_Fmfix.mat'], 'paras');

for nCell   = 1:length(totCell)  
    
    spk          = totCell(nCell).spk;
    dff          = totCell(nCell).dff;
    if isa(dff, 'single'); dff = double(dff); end
    t_frame      = totCell(nCell).CaTime;    
    fitCaTraces  = paras(nCell).fitCaTraces;

    normalized_dff     = (dff - min(dff))/(max(dff)-min(dff));
    normalized_fitCaTraces = (fitCaTraces - min(dff))/(max(dff)-min(dff));
    
    figure;
    
    subplot(3, 1, 1)    
    hold on;
    plot(t_frame, normalized_dff, '-b', 'linewid', 1);
    gridxy(spk, [], 'Color', 'k', 'linewid', 1)
    hold off
    axis off
    % xlim([t_frame(1) t_frame(end)])
    xlim([20 30])
    ylim([0 1.4])
    title('original data')

    subplot(3, 1, 2)    
    hold on;
    plot(t_frame, normalized_dff, '-b', 'linewid', 1);
    plot(t_frame, normalized_fitCaTraces, '-r', 'linewid', 1)
    gridxy(spk, [], 'Color', 'k', 'linewid', 1)
%     plot(3, 1.1, 's', 'color',  'k', 'markerfacecolor', 'k')
%     text(3.3, 1.1, 'original DF/F', 'color', 'k')
%     plot(13, 1.1, 's', 'color',  'r', 'markerfacecolor', 'r')
%     text(13.3, 1.1, 'fit', 'color', 'k')
    hold off
    axis off
    % xlim([t_frame(1) t_frame(end)])
    xlim([20 30])
    ylim([0 1.4])
    title('S2C model')
        
    setPrint(18, 10,[PlotDir 'ModelExampleCellFits/' totCell(nCell).expression '_' ...
        totCell(nCell).cellName '_' num2str(totCell(nCell).nRep,'%02d')])
    
    close all;
end

