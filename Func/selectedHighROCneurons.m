% 
% Comparison based on single unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0.2
%
% Selected the spiking cells with high ROC and firing threshold

function selectedNeuronalIndex = selectedHighROCneurons(nDataSet, params, ROCThres, selectedNeuronalIndex)
    ROCDataIndex               = filterOutLowFR(nDataSet, params, ROCThres);    
    selectedNeuronalIndex      = columnVec(selectedNeuronalIndex) & ROCDataIndex;
end


function ROCDataIndex = filterOutLowFR(nDataSet, params, ROCThres)
        
    ROCIndex           = ROCPop(nDataSet, params);
    ROCIndex           = ROCIndex > ROCThres;
    ROCDataIndex       = sum(ROCIndex(:,2:end),2)>0;
    
end