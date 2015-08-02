% 
% obtain the spike dataset from a list of files
% 
% version 2.0
%
% Comparison list
%
% Output:
% SpikeDataSet     --- yDim x 1 cells (yDims number of neurons) 
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 



function [nDataSet3D, kickOutIndexAll] = getSimultaneousCaimagingData...
                                        (nDataSet, DataSetList, ROCThres, minUnitsSession)
    
    
    ROCIndex                = DataSetList.ROCIndex;
    ActiveNeuronIndex       = DataSetList.ActiveNeuronIndex;    
    ROCDataIndex            = sum(ROCIndex(:,2:end)> ROCThres,2)>0;                                
    newDataSet              = nDataSet(ROCDataIndex & ActiveNeuronIndex);
    kickOutIndexAll         = find(~(ROCDataIndex & ActiveNeuronIndex));
    keptIndexAll            = find(ROCDataIndex & ActiveNeuronIndex);
    kickOutIndex            = [];
    nDataSet3D              = [];
    sessionIndex            = unique([newDataSet.sessionIndex]);
    unitROCIndex            = ROCIndex(keptIndexAll, :);
    
    for nSession            = 1:length(sessionIndex)
        sessionUnitIndex    = find([newDataSet.sessionIndex] == sessionIndex(nSession));
        if length(sessionUnitIndex) < minUnitsSession
            kickOutIndex    = [kickOutIndex; columnVec(sessionUnitIndex)]; %#ok<AGROW>
        else
            nDataSet3D      = [nDataSet3D, ...
                            convert2DataSet3D(newDataSet(sessionUnitIndex),...
                            unitROCIndex(sessionUnitIndex, :))]; %#ok<AGROW>
        end
    end
    
    kickOutIndexAll = [columnVec(kickOutIndexAll); ...
                        columnVec(keptIndexAll(kickOutIndex))];    
end


function nDataSet3D  = convert2DataSet3D(nDataSet, unitROCIndex)

    nDataSet3D.sessionIndex                = nDataSet(1).sessionIndex;
    nDataSet3D.nUnit                       = [nDataSet.nUnit];
    nDataSet3D.unitROCIndex                = unitROCIndex;
    nDataSet3D.unit_yes_trial_index        = nDataSet(1).unit_yes_trial_index;
    nDataSet3D.unit_no_trial_index         = nDataSet(1).unit_no_trial_index;
    
    for nUnit        = 1: length(nDataSet)
        nDataSet3D.unit_yes_trial(nUnit, :, :) = nDataSet(nUnit).unit_yes_trial;
        nDataSet3D.unit_no_trial(nUnit, :, :)  = nDataSet(nUnit).unit_no_trial;
    end

end