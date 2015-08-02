%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  generate data from shuffled dataset code (based on
%     Comparison_V1_0 && Comparison_V1_0_1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V2_1

addpath('../Func');
setDir;
load ([TempDatDir 'DataListSimultaneous.mat'], 'DataSetList');

fileToAnalysis  = [1 3 6];
ROCList         = 0.5:0.1:0.8;
minNumUnits     = 3;

for nROC                = 1:length(ROCList)
    coeffSet            = [];
    for nFile           = 1:length(fileToAnalysis)
        nData           = fileToAnalysis(nFile);
        load([TempDatDir DataSetList(nData).name '.mat']);
        
        for nSession    = 1:length(nDataSet)
            unitIndex   = nDataSet(nSession).unitROCIndex > ROCList(nROC);
            unitIndex   = sum(unitIndex(:, 2:end), 2) > 0;
            if unitIndex >= minNumUnits
                nSessionData = ...
                    [permute(nDataSet(nSession).unit_yes_trial, [2 1 3]);...
                     permute(nDataSet(nSession).unit_no_trial, [2 1 3])];
                nSessionData = normalizationDim(nSessionData, 2);
                nSessionData = nSessionData(:, unitIndex, :);
                totTargets   = ...
                    [true(length(nDataSet(nSession).unit_yes_trial_index), 1); ...
                     false(length(nDataSet(nSession).unit_no_trial_index), 1)];
                [coeffVal, gammaVal, deltaVal, correctRate] = ...
                                        coeffLDA(nSessionData, totTargets);
                
                coeffSession.fileIndex     = nData;
                coeffSession.coeffVal      = coeffVal;
                coeffSession.gammaVal      = gammaVal;
                coeffSession.deltaVal      = deltaVal;
                coeffSession.correctRate   = correctRate;
                coeffSession.sessionIndex  = nDataSet(nSession).sessionIndex;
                
                coeffSet                   = [coeffSet; coeffSession]; %#ok<AGROW>
                
            end
        end        
    end
    
    save([TempDatDir 'LDACoeffROC_0_' num2str(ROCList(nROC)*10) '.mat'], 'nDataSet');   
    
end