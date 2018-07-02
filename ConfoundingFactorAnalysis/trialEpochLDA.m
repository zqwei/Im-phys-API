%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decoding which trial period one is in, how fast can you tell the data
% that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
addNoise         = [1 0 0 0 0 0 0 0 0 0];

cmap             = cbrewer('div', 'Spectral', 128, 'cubic');

if ~exist([PlotDir '/CollectedUnitsDecodabilityEpoch'],'dir')
    mkdir([PlotDir '/CollectedUnitsDecodabilityEpoch'])
end

numRandPickUnits    = 100;
numTrials           = 500;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
trainingTargets     = trainingTargets(randperm(numTrainingTrials));
testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
testTargets         = testTargets(randperm(numTestTrials));
totTargets          = [testTargets; trainingTargets];
ROCThres            = 0.5;
numFold             = 30;

for nData           = 10 %[1 3 4]
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end
    params      = DataSetList(nData).params;
    numTime     = length(params.timeSeries);
    oldDataSet          = nDataSet;
    depth               = [DataSetList(nData).cellinfo(:).depth];
    depth               = depth(~neuronRemoveList)';
%     selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
%     selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);

    decodabilityDepth   = nan(8, numTime);

    for nDepth          = 0:100:700
        depthIndex      = depth>nDepth & depth<(nDepth + 100);        
        if sum(depthIndex)>50
            decodability    = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));    
            for nFold       = 1:numFold
                nDataSet            = oldDataSet(depthIndex);
                trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
                trainingTargets     = trainingTargets(randperm(numTrainingTrials));
                testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
                testTargets         = testTargets(randperm(numTestTrials));
                totTargets          = [testTargets; trainingTargets];
                trainingDecisions   = trainingTargets(randperm(numTrainingTrials));
                testDecisions       = testTargets(randperm(numTestTrials));
                totDecisions        = [testDecisions; trainingDecisions];
                randPickUnits       = randperm(length(nDataSet));
                randPickUnits       = randPickUnits(1:numRandPickUnits);
                nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
                decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
            end

            decodabilityDepth(nDepth/100+1, :) = mean(decodability);
        end
    end
    
    figure;
    hold on;
    imagesc(DataSetList(nData).params.timeSeries, 0:-100:-700, decodabilityDepth, [0.5 1])
    axis xy;
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    set(gca, 'yTick', -800:400:0)
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('# units');
    colormap(cmap)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'ConfoundingFactorPCA/CollectedUnitsLDADepth_' DataSetList(nData).name])
end