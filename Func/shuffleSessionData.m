%
% shuffleSessionData.m
%
%
% ----------------------------
% Output:
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function nSessionData = shuffleSessionData(nDataSet, totTargets, numTestTrials)
    
    numUnits          = length(nDataSet);
    T                 = size(nDataSet(1).unit_yes_trial, 2);
    numTrials         = length(totTargets);
    nSessionData      = zeros(numTrials, numUnits, T);
    
    for nUnit         = 1:numUnits        
        nDataSet(nUnit).unit_yes_trial = nDataSet(nUnit).unit_yes_trial(randperm(size(nDataSet(nUnit).unit_yes_trial,1)),:);
        nDataSet(nUnit).unit_no_trial  = nDataSet(nUnit).unit_no_trial(randperm(size(nDataSet(nUnit).unit_no_trial,1)),:);        
    end

    fracTest          = numTestTrials/length(totTargets);
    fracTraining      = 1 - fracTest;
    
    for nTrial        = 1:numTestTrials
        if totTargets(nTrial)
            nSessionData(nTrial, :, :)   = cell2mat(arrayfun(@(x) ...
                                          x.unit_yes_trial(...
                                          ceil(size(x.unit_yes_trial,1)*rand()*fracTest),:),...
                                          nDataSet, 'UniformOutput', false));
        else
            nSessionData(nTrial, :, :)   = cell2mat(arrayfun(@(x) ...
                                          x.unit_no_trial(...
                                          ceil(size(x.unit_no_trial,1)*rand()*fracTest),:),...
                                          nDataSet, 'UniformOutput', false));            
        end
    end
    
    for nTrial        = numTestTrials+1:numTrials
        if totTargets(nTrial)
            nSessionData(nTrial, :, :)   = cell2mat(arrayfun(@(x) ...
                                          x.unit_yes_trial(...
                                          ceil(size(x.unit_yes_trial,1)*rand()*fracTraining+...
                                          size(x.unit_yes_trial,1)*fracTest),:),...
                                          nDataSet, 'UniformOutput', false));
        else
            nSessionData(nTrial, :, :)   = cell2mat(arrayfun(@(x) ...
                                          x.unit_no_trial(...
                                          ceil(size(x.unit_no_trial,1)*rand()*fracTraining+...
                                          size(x.unit_no_trial,1)*fracTest),:),...
                                          nDataSet, 'UniformOutput', false));            
        end
    end
    
end
