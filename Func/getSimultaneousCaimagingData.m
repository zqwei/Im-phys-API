% 
% obtain the spike dataset from a list of files
% 
% version 1.1
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



function [nDataSet3D, nDataSet] = getSimultaneousCaimagingData(nDataSet, params, ROCThres, minUnitsSession)
    newDataSet             = filterOutLowFR(nDataSet, params, ROCThres);
    [nDataSet3D, nDataSet] = getSimultaneousDataSet(newDataSet, minUnitsSession);
end

    

function HighFRDataSet = filterOutLowFR(nDataSet, params, ROCThres)
        
    ROCIndex           = ROCPop(nDataSet, params);
    ROCIndex           = ROCIndex > ROCThres;
%     ROCDataIndex       = sum(ROCIndex(:,2:end),2)>2;
    ROCDataIndex       = sum(ROCIndex(:,2:end),2)>0;    
    HighFRDataSet      = nDataSet(ROCDataIndex);
    
end