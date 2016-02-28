%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Different ROC
% Same number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;

numFold             = 30;
load ([TempDatDir 'DataListShuffle.mat']);
addNoise         = [1 0 0 0];

if ~exist([PlotDir '/CollectedUnitsDecodability'],'dir')
    mkdir([PlotDir '/CollectedUnitsDecodability'])
end

cmap = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

for nData             = [1 3 4]    
    load([TempDatDir DataSetList(nData).name '.mat'])
    
    figure;
    
    %% 100 random units
    subplot(1, 4, 1);
    numRandPickUnits    = 100;
    numTrials           = numRandPickUnits*5;
    numTestTrials       = numRandPickUnits*2;
    numTrainingTrials   = numTrials - numTestTrials;
    oldDataSet          = nDataSet;
    hold on
    for ROCThres          = 0.5:0.1:0.8;
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';% & [DataSetList(nData).cellinfo(:).cellType] == 1;
        selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
        nDataSet              = oldDataSet(selectedNeuronalIndex);
        decodability          = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));        
        if numRandPickUnits> length(nDataSet); numRandPickUnits = length(nDataSet); end
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
            randPickUnits       = randPickUnits(1:numRandPickUnits);

            nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
            decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
        end
        shadedErrorBar(DataSetList(nData).params.timeSeries, mean(decodability,1),...
            std(decodability, 1)/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', cmap(ROCThres*10-4,:)}, 0.5);
    end    
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0.5 1])
    legend('boxoff')
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    set(gca, 'TickDir', 'out')
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('Decodability');
    title('100 random units')
    
    %% 200 random units
    subplot(1, 4, 2);
    
    numRandPickUnits    = 200;
    numTrials           = numRandPickUnits*5;
    numTestTrials       = numRandPickUnits*2;
    numTrainingTrials   = numTrials - numTestTrials;
    hold on
    for ROCThres          = 0.5:0.1:0.8;
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';% & [DataSetList(nData).cellinfo(:).cellType] == 1;
        selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
        nDataSet              = oldDataSet(selectedNeuronalIndex);
        decodability          = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));        
        if numRandPickUnits> length(nDataSet); numRandPickUnits = length(nDataSet); end
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
            randPickUnits       = randPickUnits(1:numRandPickUnits);

            nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
            decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
        end
        shadedErrorBar(DataSetList(nData).params.timeSeries, mean(decodability,1),...
            std(decodability, 1)/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', cmap(ROCThres*10-4,:)}, 0.5);
    end    
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0.5 1])
    legend('boxoff')
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    set(gca, 'TickDir', 'out')
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('Decodability');
    title('200 random units')
    
    %% 500 random units
    subplot(1, 4, 3);
    
    numRandPickUnits    = 500;
    numTrials           = numRandPickUnits*5;
    numTestTrials       = numRandPickUnits*2;
    numTrainingTrials   = numTrials - numTestTrials;
    hold on
    for ROCThres          = 0.5:0.1:0.8;
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';% & [DataSetList(nData).cellinfo(:).cellType] == 1;
        selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
        nDataSet              = oldDataSet(selectedNeuronalIndex);
        decodability          = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));        
        if numRandPickUnits> length(nDataSet); numRandPickUnits = length(nDataSet); end
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
            randPickUnits       = randPickUnits(1:numRandPickUnits);

            nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
            decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
        end
        shadedErrorBar(DataSetList(nData).params.timeSeries, mean(decodability,1),...
            std(decodability, 1)/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', cmap(ROCThres*10-4,:)}, 0.5);
    end    
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0.5 1])
    legend('boxoff')
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    set(gca, 'TickDir', 'out')
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('Decodability');
    title('500 random units')
    
    
    %% all units
    subplot(1, 4, 4);
    
    numRandPickUnits    = length(oldDataSet);
    numTrials           = numRandPickUnits*3;
    numTestTrials       = 500;
    numTrainingTrials   = numTrials - numTestTrials;
    hold on
    for ROCThres          = 0.5:0.1:0.8;
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';% & [DataSetList(nData).cellinfo(:).cellType] == 1;
        selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
        nDataSet              = oldDataSet(selectedNeuronalIndex);
        if numRandPickUnits> length(nDataSet); numRandPickUnits = length(nDataSet); end
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
        decodability = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
        plot(DataSetList(nData).params.timeSeries, decodability, '-', 'linewid', 1.0, 'color', cmap(ROCThres*10-4,:));
    end    
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0.5 1])
    legend('boxoff')
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    set(gca, 'TickDir', 'out')
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('Decodability');
    title('all units')
    
    %% print
    setPrint(8*4, 6, [PlotDir 'CollectedUnitsDecodability/CollectedUnitsDecodabilityROC_' DataSetList(nData).name])
end


margNames = {'0.5', '0.6', '0.7', '0.8'};

figure;
hold on
for nColor = 1:length(margNames)
    plot(0, nColor, 's', 'color', cmap(nColor,:), 'MarkerFaceColor',cmap(nColor,:),'MarkerSize', 8)
    text(1, nColor, margNames{nColor})
end
xlim([0 10])
hold off
axis off
setPrint(3, 2, [PlotDir 'CollectedUnitsDecodability/CollectedUnitsDecodabilityROC_Label'])

close all;