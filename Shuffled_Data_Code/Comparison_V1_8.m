%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_8
function Comparison_V1_8

    numTrials           = 3000;
    numTestTrials       = 600;
    numTrainingTrials   = numTrials - numTestTrials;
    trainingTargets     = rand(numTrainingTrials, 1) > 0.5;
    testTargets         = rand(numTestTrials, 1) > 0.5;
    totTargets          = [testTargets; trainingTargets];
    load ('TempDat/DataList.mat');
    addNoise            = [1 1 0 0];        
    % Comparison_V1_8_1(DataSetList, numTestTrials, totTargets, addNoise, numTrials)
    % Comparison_V1_8_2(DataSetList, numTestTrials, totTargets, addNoise, numTrials)
    Comparison_V1_8_3(DataSetList, numTestTrials, totTargets, addNoise, numTrials)
    % Comparison_V1_8_4(DataSetList, numTestTrials, totTargets, addNoise, numTrials)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.1  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Comparison_V1_8_1(DataSetList, numTestTrials, totTargets, addNoise, numTrials) %#ok<DEFNU>
    for nData             = 1:length(DataSetList)
        load(['TempDat/' DataSetList(nData).name '.mat'])        
        nSessionData = shuffleSessionData(nDataSet, totTargets);
        nSessionData = permute(nSessionData,[1 3 2]);
        EpochIndex   = epochIndex(DataSetList(nData).params);
        EpochIndex   = EpochIndex(:,ones(1,numTrials))';
        decodability = sliceDataTaper(nSessionData, EpochIndex, 0, addNoise(nData), numTestTrials);
        figure;
        hold on
        plot(DataSetList(nData).params.timeSeries, decodability, 'r', 'linewid',1);
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([0 1])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        setPrint(4, 3, ['Plot/Collected_Units_Decodability_EpochLDA1_' DataSetList(nData).name], 'png')
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.2  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched (Multi-taper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Comparison_V1_8_2(DataSetList, numTestTrials, totTargets, addNoise, numTrials) %#ok<DEFNU>
    for nData             = 1:length(DataSetList)
        load(['TempDat/' DataSetList(nData).name '.mat'])        
        nSessionData = shuffleSessionData(nDataSet, totTargets);
        nSessionData = permute(nSessionData,[1 3 2]);           
        EpochIndex   = epochIndex(DataSetList(nData).params);
        EpochIndex   = EpochIndex(:,ones(1,numTrials))';        
        figure;        
        numTaper     = 9;
        for nTaper   = 0:numTaper
            nColor       = [nTaper, 0, 0] * 0.1;
            decodability = sliceDataTaper(nSessionData, EpochIndex, nTaper, addNoise(nData), numTestTrials);
            hold on
            plot(DataSetList(nData).params.timeSeries(1:end-numTaper), decodability(1:end+nTaper-numTaper), 'Color', nColor, 'linewid',1);
        end        
        xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end-numTaper)]);
        ylim([0 1])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        setPrint(4, 3, ['Plot/Collected_Units_Decodability_EpochLDATaper1_' DataSetList(nData).name], 'png')
    end
end

function decodability               = sliceDataTaper(nSessionData, EpochIndex, taperLength, noiseData, numTestTrials)
    [numTrials, numT, numUnits]     = size(nSessionData);
    nSessionData                    = nSessionData(:,taperLength+1:end,:);
    EpochIndex                      = EpochIndex(:,1:end-taperLength);
    numT                            = numT - taperLength;
    nSessionData                    = reshape(nSessionData, numTrials*numT, numUnits); 
    EpochIndex                      = reshape(EpochIndex, numTrials*numT, 1);
    testEpoch                       = EpochIndex(1:numTestTrials*numT);
    trainingEpoch                   = EpochIndex(numTestTrials*numT+1:end);
    decodability                    = decodeEpochLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* noiseData, trainingEpoch, testEpoch, numT);    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.3  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched (Multi-LDA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Comparison_V1_8_3(DataSetList, numTestTrials, totTargets, addNoise, numTrials) %#ok<DEFNU>
    for nData             = 1:length(DataSetList)
        load(['TempDat/' DataSetList(nData).name '.mat'])        
        nSessionData = shuffleSessionData(nDataSet, totTargets);
        nSessionData = permute(nSessionData,[1 3 2]);           
        EpochIndex   = epochIndex(DataSetList(nData).params);
        EpochIndex   = EpochIndex(:,ones(1,numTrials))';        
        timePoints   = timePointTrialPeriod(DataSetList(nData).params.polein, DataSetList(nData).params.poleout, DataSetList(nData).params.timeSeries);    
        figure;
        hold on;    
        for nPeriods        = 1: length(timePoints) -2 
            nColor          = [nPeriods, 0, 0] * 0.3;
            sliceIndex      = timePoints(nPeriods):timePoints(nPeriods+2);
            decodability    = sliceData(nSessionData, EpochIndex, sliceIndex, addNoise(nData), numTestTrials);
            plot(DataSetList(nData).params.timeSeries(sliceIndex), decodability, 'Color', nColor, 'linewid',1);
        end
        xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
        ylim([0 1])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color', 'k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        setPrint(4, 3, ['Plot/Collected_Units_Decodability_EpochLDA3_' DataSetList(nData).name], 'png')
    end
end

function decodability               = sliceData(nSessionData, EpochIndex, sliceIndex, noiseData, numTestTrials)
    [numTrials, ~, numUnits]        = size(nSessionData);
    nSessionData                    = nSessionData(:,sliceIndex,:);
    EpochIndex                      = EpochIndex(:,sliceIndex);
    numT                            = length(sliceIndex);
    nSessionData                    = reshape(nSessionData, numTrials*numT, numUnits); 
    EpochIndex                      = reshape(EpochIndex, numTrials*numT, 1);
    testEpoch                       = EpochIndex(1:numTestTrials*numT);
    trainingEpoch                   = EpochIndex(numTestTrials*numT+1:end);
    decodability                    = decodeEpochLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* noiseData, trainingEpoch, testEpoch, numT);    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.4  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched (Multi-LDA Multi-Taper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
