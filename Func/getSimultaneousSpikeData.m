% 
% obtain the spike dataset from a list of files
% 
% version 2.0
%
% Comparison list
%
% Based on plotSimultaneousSpikeData
%
% Output:
% SpikeDataSet     --- yDim x 1 cells (yDims number of neurons) 
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 



function [nDataSet3D, kickOutIndexAll] = getSimultaneousSpikeData(...
                                    nDataSet, DataSetList, minRate, perMinRate, ...
                                    ROCThres, minUnitsSession)
    
    ROCIndex          = DataSetList.ROCIndex;
    ActiveNeuronIndex = DataSetList.ActiveNeuronIndex;                        
                                
    %%%
    % Set minimal ratio between number of trial and number of unit
    %%%
    trialUnitRatio     = 2.5;
    minNumTrial        = 20;
    
    %%%
    % Find high active units, where the probability of spike is > perMinRate
    %%%
    HighFRDataSetIndex = arrayfun(@(x) mean(mean([x.unit_yes_trial; x.unit_no_trial])>minRate)>perMinRate, ...
                         nDataSet, 'Uniformoutput', false);                                
    HighFRDataSetIndex = cell2mat(HighFRDataSetIndex);
    
    %%%
    % Find high selective units, where the probability of ROC is > ROCThres
    %%%
    ROCDataIndex       = sum(ROCIndex(:,2:end) > ROCThres,2)>0;
    kickOutIndexAll    = find(~(ActiveNeuronIndex & HighFRDataSetIndex & ROCDataIndex));
    keptIndexAll       = find(ActiveNeuronIndex & HighFRDataSetIndex & ROCDataIndex);
    newDataSet         = nDataSet(ActiveNeuronIndex & HighFRDataSetIndex & ROCDataIndex);
    unitROCIndex       = ROCIndex(keptIndexAll, :);
    
    sessionIndex            = unique([newDataSet.sessionIndex]);
    nDataSet3D              = [];
    kickOutIndex            = [];
    
    for nSession            = 1:length(sessionIndex)
        sessionUnitIndex    = find([newDataSet.sessionIndex] == sessionIndex(nSession));
        
%         if sessionIndex(nSession)  == 75; keyboard; end
        
        if length(sessionUnitIndex) < minUnitsSession
            kickOutIndex    = [kickOutIndex; columnVec(sessionUnitIndex)]; %#ok<AGROW>
        else
            trialUnitMat    = zeros(32, 600);
            for nUnit       = 1:length(sessionUnitIndex)
                trialUnitMat(newDataSet(sessionUnitIndex(nUnit)).nUnit, ...
                            newDataSet(sessionUnitIndex(nUnit)).unit_yes_trial_index)...
                            = 1;
                trialUnitMat(newDataSet(sessionUnitIndex(nUnit)).nUnit, ...
                            newDataSet(sessionUnitIndex(nUnit)).unit_no_trial_index)...
                            = -1;
            end
            trialIndex      = sum(trialUnitMat~=0, 1) >= minUnitsSession;
            unitIndex       = sum(trialUnitMat(:, trialIndex)==1, 2)> minNumTrial...
                            & sum(trialUnitMat(:, trialIndex)==-1, 2)> minNumTrial;
            subUnitIndex    = ismember([newDataSet(sessionUnitIndex).nUnit], find(unitIndex));
            kickOutIndex    = [kickOutIndex; columnVec(sessionUnitIndex(~subUnitIndex))]; %#ok<AGROW>
            sessionUnitIndex = sessionUnitIndex(subUnitIndex);
            if sum(unitIndex) < minUnitsSession
                kickOutIndex   = [kickOutIndex; columnVec(sessionUnitIndex)]; %#ok<AGROW>
            else
                
                trialIndex     = sum(trialUnitMat(unitIndex, :)~=0, 1) == sum(unitIndex);
                
                if sum(trialIndex & sum(trialUnitMat,1)<0) > sum(unitIndex) * trialUnitRatio ...
                    && sum(trialIndex & sum(trialUnitMat,1)>0) > sum(unitIndex) * trialUnitRatio
                    % Keep the whole section in DataSet3D
                    nDataSet3D  = [nDataSet3D, ...
                                convert2DataSet3D(newDataSet(sessionUnitIndex),...
                                unitROCIndex(sessionUnitIndex, :), find(trialIndex))]; %#ok<AGROW>
                else % keep one unit out
                    newTrialUnitMat     = trialUnitMat(unitIndex, :);
                    kickUnit            = 0;
                    numTrialAfterKick   = 0;
                    trialIndexAfterKick = false(1, size(newTrialUnitMat, 2));
                    for nKickUnit       = 1:size(newTrialUnitMat)
                        subTrialUnitMat = newTrialUnitMat;
                        subTrialUnitMat(nKickUnit, :) = [];
                        trialIndex      = sum(subTrialUnitMat~=0, 1) == size(subTrialUnitMat, 1);
                        
                        if numTrialAfterKick < sum(trialIndex)
                            numTrialAfterKick   = sum(trialIndex);
                            kickUnit            = nKickUnit;
                            trialIndexAfterKick = trialIndex;
                        end
                    end
                    
                    if kickUnit         == 0
                        kickOutIndex   = [kickOutIndex; columnVec(sessionUnitIndex)]; %#ok<AGROW>
                    else
                        subTrialUnitMat = newTrialUnitMat;
                        subTrialUnitMat(kickUnit, :) = [];
                        kickOutUnitIndex = sessionUnitIndex;
                        kickOutUnitIndex(kickUnit)   = [];
                        
                        if sum(trialIndexAfterKick & sum(subTrialUnitMat,1)<0) > size(subTrialUnitMat, 1) * trialUnitRatio ...
                            && sum(trialIndexAfterKick & sum(subTrialUnitMat,1)>0) > size(subTrialUnitMat, 1) * trialUnitRatio
                            kickOutIndex   = [kickOutIndex; sessionUnitIndex(kickUnit)]; %#ok<AGROW>
                            % Keep the kick-one-out section in DataSet3D
                            nDataSet3D  = [nDataSet3D, ...
                                        convert2DataSet3D(newDataSet(kickOutUnitIndex),...
                                        unitROCIndex(kickOutUnitIndex, :), find(trialIndexAfterKick))]; %#ok<AGROW>
                        else
                            kickOutIndex   = [kickOutIndex; columnVec(sessionUnitIndex)]; %#ok<AGROW>
                        end
                    end
                    
                end
            end
        end
    end
    
    kickOutIndexAll = [columnVec(kickOutIndexAll); ...
                        columnVec(keptIndexAll(kickOutIndex))];
    
end


function nDataSet3D  = convert2DataSet3D(nDataSet, unitROCIndex, trialIndex)

%     disp(nDataSet(1).sessionIndex)

    nDataSet3D.sessionIndex                = nDataSet(1).sessionIndex;
    nDataSet3D.nUnit                       = [nDataSet.nUnit];
    nDataSet3D.unitROCIndex                = unitROCIndex;
    nDataSet3D.unit_yes_trial_index        = nDataSet(1).unit_yes_trial_index...
                                            (ismember(nDataSet(1).unit_yes_trial_index, trialIndex));
    nDataSet3D.unit_no_trial_index         = nDataSet(1).unit_no_trial_index...
                                            (ismember(nDataSet(1).unit_no_trial_index, trialIndex));
    
    for nUnit        = 1: length(nDataSet)
        nDataSet3D.unit_yes_trial(nUnit, :, :) = nDataSet(nUnit).unit_yes_trial...
                                                (ismember(nDataSet(nUnit).unit_yes_trial_index, trialIndex), :);
        nDataSet3D.unit_no_trial(nUnit, :, :)  = nDataSet(nUnit).unit_no_trial...
                                                (ismember(nDataSet(nUnit).unit_no_trial_index, trialIndex), :);
    end
    
    nDataSet3D.unit_yes_trial               = permute(nDataSet3D.unit_yes_trial, [2 1 3]);
    nDataSet3D.unit_no_trial                = permute(nDataSet3D.unit_no_trial, [2 1 3]);
end