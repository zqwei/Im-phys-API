function firingRates        = generateDPCADataVec(nDataSet, numTrials)

    numUnits          = length(nDataSet);
    T                 = size(nDataSet(1).unit_yes_trial, 2);
    numStim           = 2;
    firingRates       = zeros(numUnits, numStim, T, numTrials);
    
    for nTrial        = 1:numTrials
        firingRates(:, 1, :, nTrial)  = cell2mat(arrayfun(@(x) ...
                                      x.unit_yes_trial(...
                                      ceil(size(x.unit_yes_trial,1)*rand()),:)',...
                                      nDataSet, 'UniformOutput', false))';
        firingRates(:, 2, :, nTrial)  = cell2mat(arrayfun(@(x) ...
                                      x.unit_no_trial(...
                                      ceil(size(x.unit_no_trial,1)*rand()),:)',...
                                      nDataSet, 'UniformOutput', false))';            
    end


end
