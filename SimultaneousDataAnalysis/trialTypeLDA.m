%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Same ROC
% Different number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('../Func');
setDir;

load ([TempDatDir 'DataListSimultaneous.mat']);

if ~exist([PlotDir '/SimultaneousUnitsDecodability'],'dir')
    mkdir([PlotDir '/SimultaneousUnitsDecodability'])
end



for nData             = 1:length(DataSetList)
    load([TempDatDir 'SimultaneousLDA_' DataSetList(nData).name '.mat'], 'correctRates');
    load([TempDatDir 'SimultaneousUnitStats_' DataSetList(nData).name '.mat'], 'numSimultaneousUnits', 'numSimultaneousTrial');
    
    [ss, sessionIndex] = sortrows([numSimultaneousUnits, numSimultaneousTrial], [-1 -2]);
    

    
    
    figure;
    
    for nSession      = 1:length(correctRates)
        correctRate   = correctRates{sessionIndex(nSession)};        
        subplot(m, m, nSession);
        plot(DataSetList(nData).params.timeSeries, correctRate, 'Color','k','Linestyle','-','linewid', 2.0);
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([0.5 1])
        legend('boxoff')
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        xlabel('Time (s)');
        ylabel('Decodability');
        title(['# units: ' num2str(ss(nSession, 1)) '; # trials: ' num2str(ss(nSession, 2))])
    end
    
    
    
    
    
    
    setPrint(8*m, 6*m, [PlotDir 'SimultaneousUnitsDecodability/SimultaneousUnitsDecodability_' DataSetList(nData).name], 'pdf')
end


close all;









numFold             = 10;
numTrials           = 1000;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
ROCThres            = 0.5;

load ([TempDatDir 'DataListShuffle.mat']);
addNoise            = [1 0 0 0];

cmap                = cbrewer('div', 'Spectral', 128, 'cubic');


if ~exist([PlotDir '/CollectedUnitsDecodability'],'dir')
    mkdir([PlotDir '/CollectedUnitsDecodability'])
end

stepSize            = 10;

for nData             = [1 3 4]
    load([TempDatDir DataSetList(nData).name '.mat'])
    oldDataSet               = nDataSet;
    maxRandPickUnits         = 20;
    decodabilityAll          = zeros(maxRandPickUnits, size(nDataSet(1).unit_yes_trial,2));
    figure;
    hold on
    for numRandPickUnits      = 1:maxRandPickUnits;
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';% & [DataSetList(nData).cellinfo(:).cellType] == 1;
        selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
        nDataSet              = oldDataSet(selectedNeuronalIndex);
        decodability          = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));

        for nFold    = 1:numFold
            trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
            trainingTargets     = trainingTargets(randperm(numTrainingTrials));
            testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
            testTargets         = testTargets(randperm(numTestTrials));
            totTargets          = [testTargets; trainingTargets];

            trainingDecisions   = trainingTargets(randperm(numTrainingTrials));
            testDecisions       = testTargets(randperm(numTestTrials));
            totDecisions        = [testDecisions; trainingDecisions];

            randPickUnits       = randperm(length(nDataSet));
            randPickUnits       = randPickUnits(1:numRandPickUnits*stepSize);

            nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
            decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
        end
        decodabilityAll(numRandPickUnits, :) = mean(decodability,1);
    end
    imagesc(DataSetList(nData).params.timeSeries, (1:maxRandPickUnits)*stepSize, decodabilityAll,[0.5 1]);
    axis xy;
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([stepSize maxRandPickUnits*stepSize])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('# units');
    colormap(cmap)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'CollectedUnitsDecodability/CollectedUnitsDecodabilityFixedROCThres_0_5_' DataSetList(nData).name])
end


ROCThres            = 0.7;

for nData             = [1 3 4]
    load([TempDatDir DataSetList(nData).name '.mat'])
    oldDataSet               = nDataSet;
    maxRandPickUnits         = 20;
    decodabilityAll          = zeros(maxRandPickUnits, size(nDataSet(1).unit_yes_trial,2));
    figure;
    hold on
    for numRandPickUnits      = 1:maxRandPickUnits;
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';% & [DataSetList(nData).cellinfo(:).cellType] == 1;
        selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
        nDataSet              = oldDataSet(selectedNeuronalIndex);
        decodability          = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));

        for nFold    = 1:numFold
            trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
            trainingTargets     = trainingTargets(randperm(numTrainingTrials));
            testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
            testTargets         = testTargets(randperm(numTestTrials));
            totTargets          = [testTargets; trainingTargets];

            trainingDecisions   = trainingTargets(randperm(numTrainingTrials));
            testDecisions       = testTargets(randperm(numTestTrials));
            totDecisions        = [testDecisions; trainingDecisions];

            randPickUnits       = randperm(length(nDataSet));
            randPickUnits       = randPickUnits(1:numRandPickUnits*stepSize);

            nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
            decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
        end
        decodabilityAll(numRandPickUnits, :) = mean(decodability,1);
    end
    imagesc(DataSetList(nData).params.timeSeries, (1:maxRandPickUnits)*stepSize, decodabilityAll,[0.5 1]);
    axis xy;
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([stepSize maxRandPickUnits*stepSize])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('# units');
    colormap(cmap)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'CollectedUnitsDecodability/CollectedUnitsDecodabilityFixedROCThres_0_7_' DataSetList(nData).name])
end

setColorbar(cmap, 0.5, 1, 'accuracy', [PlotDir 'CollectedUnitsDecodability/CollectedUnitsDecodabilityFixedROCThres_0_7_'])
setColorbar(cmap, 0.5, 1, 'accuracy', [PlotDir 'CollectedUnitsDecodability/CollectedUnitsDecodabilityFixedROCThres_0_5_'])


close all;