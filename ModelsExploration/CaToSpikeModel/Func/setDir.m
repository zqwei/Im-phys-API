% setDir.m
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


% function setDir

% dataDir     = '../../../Data_In_Use/Forward_Model/';
% dataFolders = {'GCaMP6f_11cells (Chen 2013)/', 'GCaMP6s_9cells (Chen 2013)/'};
% expression  = {'virus', 'virus'};
% CaIndicator = {'GCaMP6f', 'GCaMP6s'};

TempDir     = '../TempDat/';
if ~exist(TempDir,'dir')
    mkdir(TempDir)
end

PlotDir     = '../Plots/';
if ~exist(PlotDir,'dir')
    mkdir(PlotDir)
end