%
% extractALMImagingData.m
%
% This performs the transformation from the raw ROI data to the Ca++ for
% KD's data
% 
% imaging data (dFF)
% 
% Input:
% ROI_list     --- structure of raw Ca++ imaging data
%
% ----------------------------
% Output:
% dFF    --- yDim x T x K tensor (3D) of data
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function [dFF, mean_F, correctRightTrial, correctLeftTrial, errorRightTrial, errorLeftTrial, cell_position] = ...
    extractALMImagingDataKD (dat, params)
    % dat is dat_small in KD's file
    
    ROI_list = dat.roi;
    
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    cueTime     = 3.13;

    neuronNum          = length(ROI_list);
    [numTime, trialNum]= size(ROI_list(1).intensity);
    
    if numTime > 135
        keyboard()
    end
    
    neuronSig          = nan(neuronNum, floor(numTime/2), trialNum);
    neuropilSig        = nan(neuronNum, floor(numTime/2), trialNum);
    
    cell_position = reshape([ROI_list.centerXY], 2, neuronNum)';
    
    for nNeuron        = 1:neuronNum
        neuronSig(nNeuron, :, :)   = ROI_list(nNeuron).intensity(1:2:133, :)+ROI_list(nNeuron).intensity(2:2:134, :);
        neuropilSig(nNeuron, :, :) = ROI_list(nNeuron).neuropil_intensity(1:2:133, :)+ROI_list(nNeuron).neuropil_intensity(2:2:134, :);
    end
    
    correctedNeuronSig = neuronSig - 0.7*neuropilSig;
    correctedNeuronSig = correctedNeuronSig/2;
    
    correct           = dat.trial_type == dat.lick_direction;
    type              = dat.trial_type;
    early             = dat.early_lick==1;
    if length(early) > length(correct)
        early         = early(1:length(correct));
    end
    correctRightTrial = correct & (type==1) & (~early);
    correctLeftTrial  = correct & (type==2) & (~early);
    errorRightTrial   = ~correct & (type==1) & (~early);
    errorLeftTrial    = ~correct & (type==2) & (~early);
    
    alignedNeuralData = correctedNeuronSig(:, round((cueTime+params.polein-0.5)*params.frameRate):round((cueTime+1.2)*params.frameRate),:);
    timePointNum      = timePoints(end);
    
    mean_F             = squeeze(mean(mean(alignedNeuralData(:,timePoints(1):timePoints(2),:),3),2));
    alignedNeuralData  = alignedNeuralData - repmat(mean_F,1,timePointNum,trialNum);
    dFF                = alignedNeuralData./(repmat(mean_F,1,timePointNum,trialNum)+0.0001);
