% 
% Comparison based on single unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0.2
%
% Selected the spiking cells with high ROC and firing threshold

function Comparison_V1_0_2
    addpath('../Func');
    setDir;
    load ([TempDatDir 'DataListModeled.mat']);

    % -------------------------------------------------------------------------
    % Spiking Dataset
    % -------------------------------------------------------------------------

    nData                       = 1;
    load([TempDatDir DataSetList(nData).name '.mat']);
    % only considering pyramidal cells
    selectedNeuronalIndex       = [DataSetList(nData).cellinfo(:).cellType]==1;
    minRate                     = 5; % Hz
    perMinRate                  = 0.4;
    ROCThres                    = 0.6;
    HighFRDataSet               = filterOutLowFR(nDataSet, DataSetList(nData).params, minRate, perMinRate, ROCThres);
    selectedNeuronalIndex       = selectedNeuronalIndex' & HighFRDataSet; %#ok<NASGU>
    save([TempDatDir 'DataListModeled.mat'], 'selectedNeuronalIndex', '-append');
end


function HighFRDataSet = filterOutLowFR(nDataSet, params, minRate, perMinRate, ROCThres)
    
    HighFRDataSetIndex = arrayfun(@(x) mean(mean([x.unit_yes_trial; x.unit_no_trial])>minRate)>perMinRate, ...
                         nDataSet, 'Uniformoutput', false);
    HighFRDataSetIndex = cell2mat(HighFRDataSetIndex);
    
    ROCIndex           = ROCPop(nDataSet, params);
    ROCIndex           = ROCIndex > ROCThres;
    ROCDataIndex       = sum(ROCIndex(:,2:end),2)>0;
    
    HighFRDataSet      = HighFRDataSetIndex & ROCDataIndex;
    
end