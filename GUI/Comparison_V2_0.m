% 
% Comparison based on single unit acitivity (GUI-version)
%
% -------------------------------------------------------------------------
% version 1.0
% + All datasets are loaded at the same time
% 
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

% Load all datasets
% Comparison_V1_0

% Check TempDat Directory
% TempDat Directory -- most of the initialized information

if ~exist('TempDat','dir')
    mkdir('TempDat')
end

% Load GUI
Load_DataSets