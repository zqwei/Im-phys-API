%
% shuffleSessionDataWithError.m
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

function nSessionData = shuffleSessionDataWithError(nDataSet, totTargets, totDecision)
    
    numUnits          = length(nDataSet);
    T                 = size(nDataSet(1).unit_yes_trial, 2);
    numTrials         = length(totTargets);
    nSessionData      = zeros(numTrials, numUnits, T);
    
    
    for nTrial        = 1:numTrials
        if totTargets(nTrial) && totTargets(nTrial) == totDecision(nTrial)
            nSessionData(nTrial, :, :)   = cell2mat(arrayfun(@(x) ...
                                          x.unit_yes_trial(...
                                          ceil(size(x.unit_yes_trial,1)*rand()),:),...
                                          nDataSet, 'UniformOutput', false));
        elseif ~totTargets(nTrial) && totTargets(nTrial) == totDecision(nTrial)
            nSessionData(nTrial, :, :)   = cell2mat(arrayfun(@(x) ...
                                          x.unit_no_trial(...
                                          ceil(size(x.unit_no_trial,1)*rand()),:),...
                                          nDataSet, 'UniformOutput', false));   
        elseif totTargets(nTrial) && totTargets(nTrial) ~= totDecision(nTrial)
            nSessionData(nTrial, :, :)   = cell2mat(arrayfun(@(x) ...
                                          x.unit_yes_error(...
                                          ceil(size(x.unit_yes_error,1)*rand()),:),...
                                          nDataSet, 'UniformOutput', false));
        elseif ~totTargets(nTrial) && totTargets(nTrial) ~= totDecision(nTrial)
            nSessionData(nTrial, :, :)   = cell2mat(arrayfun(@(x) ...
                                          x.unit_no_error(...
                                          ceil(size(x.unit_no_error,1)*rand()),:),...
                                          nDataSet, 'UniformOutput', false));
        end
    end
    
end
