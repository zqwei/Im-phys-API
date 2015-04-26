%
% extractALMImagingData.m
%
% This performs the transformation from the raw ROI data to the Ca++
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


function [dFF, correctRightTrial, correctLeftTrial, errorRightTrial, errorLeftTrial] = ...
    extractALMImagingDataRaw (ROI_list, trial, nimage, params)
    
    if isfield(params, 'timeWindowIndexRange')
        timeWindowIndexRange = params.timeWindowIndexRange;  % sec
    else
        timeWindowIndexRange = -45:30;  % sec
    end    

    neuronNum          = length(ROI_list);
    trialNum           = length(trial);
    timePointNum       = length(timeWindowIndexRange);


    neuronSig          = [ROI_list(:).fmean];
    correctedNeuronSig = neuronSig;
    validTrials        = nimage > 0;

    cueFrame           = round(params.frameRate*[trial.cuetime]);
    trialIndexOffset   = [0,cumsum(nimage(1:(end-1)))];
    cueFrame           = cueFrame+trialIndexOffset;

    % dff=(fmean-f0)/f0;
    % falign=zeros(length(xrange),trialNum);
    alignedNeuralData  = zeros(neuronNum,timePointNum,trialNum);

    for ii=1:length(cueFrame)
      alignedNeuralData(:,:,ii) = correctedNeuronSig(cueFrame(ii)+timeWindowIndexRange,:)';
    end

    % presampleFrame      = cueFrame'*ones(1,length(presample)) + ones(length(cueFrame),1)*presample;
    
    % mean_F              = mean(correctedNeuronSig(presampleFrame(:),:),1)';
    % mean_F              = mean(correctedNeuronSig)';
    
    % mean_F             = squeeze(mean(mean(alignedNeuralData,3),2));
    
    correct           = [trial.correct];
    type              = [trial.type];
    early             = [trial.early];
    correctRightTrial = correct & (type=='r') & (~early) & validTrials;
    correctLeftTrial  = correct & (type=='l') & (~early) & validTrials;
    errorRightTrial   = ~correct & (type=='r') & (~early) & validTrials;
    errorLeftTrial    = ~correct & (type=='l') & (~early) & validTrials;
    
    dFF                = alignedNeuralData;

