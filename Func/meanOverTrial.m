%
% meanOverTrial.m
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

function tmeanData            = meanOverTrial(nSessionData)
    if length(nSessionData)   == 1
        nUnits                = length(nSessionData.nUnit);
        numT                  = size(nSessionData.unit_yes_trial, 3);
%         meanData              = zeros(nUnits, 2, 1, numT);
        tmeanData             = zeros(nUnits, 2, numT);
        tmeanData(:,1, :)     = mean(nSessionData.unit_yes_trial, 2);
        tmeanData(:,2, :)     = mean(nSessionData.unit_no_trial, 2);
%         meanData(:,:,1,:)     = tmeanData;
    else
        nUnits                = length(nSessionData);
        numT                  = size(nSessionData(1).unit_yes_trial, 2);
%         meanData              = zeros(nUnits, 2, 1, numT);
        tmeanData             = zeros(nUnits, 2, numT);
        unitMean              = arrayfun(@(x) mean(x.unit_yes_trial), nSessionData, 'UniformOutput',false);
        tmeanData(:,1, :)     = cell2mat(unitMean);
        unitMean              = arrayfun(@(x) mean(x.unit_no_trial), nSessionData, 'UniformOutput',false);
        tmeanData(:,2, :)     = cell2mat(unitMean);
%         meanData(:,:,1,:)     = tmeanData;
    end
end