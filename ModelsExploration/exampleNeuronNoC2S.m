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

load([TempDatDir 'ParamsFitCells_S2CModel_Fmfix.mat'], 'paras');

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
        
%     figure;
%     hold on;
%     plot(t_frame, normalized_dff, '-k', 'linewid', 1);
%     plot(spk, ones(size(spk))*1.05, '+k', 'linewid', 1)
%     plot(t_frame, normalized_fitCaTraces, '--', 'color',  [0.6350    0.0780    0.1840], 'linewid', 1)
%     hold off
%     axis off
%     xlim([t_frame(1) t_frame(end)])
% %     xlim([0 30])
%     ylim([0 1.4])
%     title('original data')
% 
%     setPrint(18, 10,[PlotDir 'ModelExampleCellFits/' totCell(nCell).expression '_' ...
%         totCell(nCell).cellName '_' num2str(totCell(nCell).nRep,'%02d_noC2S')])
%     setPrint(18, 10,[PlotDir 'ModelExampleCellFits/' totCell(nCell).expression '_' ...
%         totCell(nCell).cellName '_' num2str(totCell(nCell).nRep,'%02d_noC2S')], 'png')
    
    figure;
    [fout, xout] = ksdensity(normalized_dff - normalized_fitCaTraces - mean(normalized_dff - normalized_fitCaTraces));
    [~, pValue ] = kstest(normalized_dff - normalized_fitCaTraces - mean(normalized_dff - normalized_fitCaTraces));
    plot(xout, fout, '-', 'linewid', 1);
    box off
    xlabel('dff - dff(fit)')
    ylabel('prob. dens.')
    title(['ks-test p = ' num2str(pValue, '%.3f')])

    setPrint(8, 6,[PlotDir 'ModelExampleCellFits/' totCell(nCell).expression '_' ...
        totCell(nCell).cellName '_' num2str(totCell(nCell).nRep,'%02d_noC2S_Distr')])
    setPrint(8, 6,[PlotDir 'ModelExampleCellFits/' totCell(nCell).expression '_' ...
        totCell(nCell).cellName '_' num2str(totCell(nCell).nRep,'%02d_noC2S_Distr')], 'png')
    close all;
    
end

