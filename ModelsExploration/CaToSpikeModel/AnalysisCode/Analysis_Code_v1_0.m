%
% Comparison between linear model and TWC model
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

% function Analysis_Code_v1_0 

addpath('../Func/')
setDir;
load([TempDir 'TotCell.mat'], 'totCell');

for tCell         = 1:length(totCell)
    para.t_frame                 = totCell(tCell).CaTime;
    para.t_ephys                 = totCell(tCell).ephysTime;
    para.fneuropil               = totCell(tCell).fNeuropil;
    para.fmean                   = totCell(tCell).fROI;
    para.filt                    = totCell(tCell).filteredEphys;
    para.peak                    = totCell(tCell).detectedSpikes;
    dff                          = get_baseline(para);
    normalized_dff               = (dff - min(dff))/(max(dff) - min(dff));
    SAMP                         = cont_ca_sampler(normalized_dff);
    plot_continuous_samples(SAMP, normalized_dff, para);
end