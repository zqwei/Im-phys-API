%
% Collecting parameters from TCW's file
% Includes:
% 1. Fm
% 2. K
% 3. n
% 4. tau_rise
% 5. tau_decay
%
% Plot ---
%
% ca++ trace     -b
% spike time     +k
% modeled trace  -r
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%  

addpath('../Func/')
setDir;

load([TempDir 'TotCell.mat'], 'totCell');

load([TempDir 'Paras_TWC_model.mat'], 'paras','cellToAnalysis'); 

if ~exist([PlotDir 'TWCModelFits/'],'dir')
    mkdir([PlotDir 'TWCModelFits/'])
end

numCell                             = 0;


for nCell   = 1:length(totCell)
    if cellToAnalysis(nCell)
        numCell                      = numCell + 1;
        para.t_frame                 = totCell(nCell).CaTime;
        para.t_ephys                 = totCell(nCell).ephysTime;
        para.fneuropil               = totCell(nCell).fNeuropil;
        para.fmean                   = totCell(nCell).fROI;
        para.filt                    = totCell(nCell).filteredEphys;
        para.peak                    = totCell(nCell).detectedSpikes;
        dff                          = get_baseline(para);
        spk                          = {para.t_ephys(para.peak)};
        para_final                   = [paras(numCell).Fm paras(numCell).K paras(numCell).n paras(numCell).tau_d paras(numCell).tau_r];
        CaTraces                     = func_getCaTraces_general_new(spk,para.t_frame,para_final);
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
      