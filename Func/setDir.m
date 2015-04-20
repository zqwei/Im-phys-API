%
% setDir.m
%
% This file sets up the basic direction information of the Ca++ imaging
% data and the spiking data
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

warning('off', 'all')
% addpath('./Func')

% function setDir

SpikingDataDir              = '../../../Data_In_Use/Dataset_Comparison/ElectrophysiologyData/';
CaImagingDataDir            = '../../../Data_In_Use/Dataset_Comparison/ImagingData/';
CaImagingShortDelayFastDir  = [CaImagingDataDir 'delay1f4t/'];
CaImagingShortDelaySlowDir  = [CaImagingDataDir 'delay1s4t/'];
CaImagingShortDelaySlowVirusDir  = [CaImagingDataDir 'delay1s4v/'];
CaImagingLongDelayFastDir   = [CaImagingDataDir 'delay3f0t/'];
CaImagingLongDelaySlowDir   = [CaImagingDataDir 'delay3s0t/'];


% Load filenames
SpikeFileList                   = dir([SpikingDataDir '*.mat']);
CaImagingShortDelayFastFileList = dir([CaImagingDataDir 'delay1f4t/*.mat']);
CaImagingShortDelaySlowFileList = dir([CaImagingDataDir 'delay1s4t/*.mat']);
CaImagingShortDelaySlowVirusFileList  = dir([CaImagingDataDir 'delay1s4v/*.mat']);
CaImagingLongDelayFastFileList  = dir([CaImagingDataDir 'delay3f0t/*.mat']);
CaImagingLongDelaySlowFileList  = dir([CaImagingDataDir 'delay3s0t/*.mat']);

TempDatDir                   = '../TempDat/';
PlotDir                      = '../Plot/';