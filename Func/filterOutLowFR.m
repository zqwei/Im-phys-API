function HighFRDataSet = filterOutLowFR(nDataSet, params, minRate, perMinRate, ROCThres)
    
    HighFRDataSetIndex = arrayfun(@(x) mean(mean([x.unit_yes_trial; x.unit_no_trial])>minRate)>perMinRate, ...
                         nDataSet, 'Uniformoutput', false);
    HighFRDataSetIndex = cell2mat(HighFRDataSetIndex);
    
    ROCIndex           = ROCPop(nDataSet, params);
    ROCIndex           = ROCIndex > ROCThres;
    ROCDataIndex       = sum(ROCIndex(:,2:end),2)>0;
    
    HighFRDataSet      = HighFRDataSetIndex & ROCDataIndex;
    
end