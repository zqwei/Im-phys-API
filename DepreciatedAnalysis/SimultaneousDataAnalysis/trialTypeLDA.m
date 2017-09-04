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

mCol             = 4;
addNoise         = [1 0 0 0];

for nData        = [1 3 4]
    
    %% get #s of trials and units in each datasets    
    load([TempDatDir DataSetList(nData).name '.mat']);
    mRow = ceil(length(nDataSet)/mCol);
    
    figure;
    
    for nSession  = 1:length(nDataSet)
        
        %% compute decodability of simultaneous recording data
        
        numYesTrial = length(nDataSet(nSession).unit_yes_trial_index);
        numNoTrial  = length(nDataSet(nSession).unit_no_trial_index);
        totTargets  = [true(numYesTrial, 1); false(numNoTrial, 1)];
        numFold     = 10;
        numUnits   = length(nDataSet(nSession).nUnit);
        numTrials  = numYesTrial + numNoTrial;
        
        simDecodability = zeros(numFold, size(nDataSet(nSession).unit_yes_trial, 3));
        for nFold   = 1:numFold
            testIndex = [randperm(numYesTrial, 4)'; randperm(numNoTrial, 4)'+numYesTrial];
            trainingIndex    = true(size(totTargets));
            trainingIndex(testIndex) = false;
            trainingIndex    = find(trainingIndex);
            trainingTargets  = totTargets(trainingIndex);
            testTargets      = totTargets(testIndex);
            nSessionData     = [nDataSet(nSession).unit_yes_trial; nDataSet(nSession).unit_no_trial];
            nSessionData     = nSessionData([testIndex; trainingIndex], :, :);
            simDecodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
        end
        
        %% compute decodability of shuffled recording data
        shfDecodability = zeros(numFold, size(nDataSet(nSession).unit_yes_trial, 3));
        for nFold   = 1:numFold
            testIndex = [randperm(numYesTrial, 4)'; randperm(numNoTrial, 4)'+numYesTrial];
            trainingIndex    = true(size(totTargets));
            trainingIndex(testIndex) = false;
            trainingIndex    = find(trainingIndex);
            trainingTargets  = totTargets(trainingIndex);
            testTargets      = totTargets(testIndex);
            nSessionData     = [nDataSet(nSession).unit_yes_trial; nDataSet(nSession).unit_no_trial];  
            nSessionTrainData = nSessionData(trainingIndex, :, :);
            nSessionTestData  = nSessionData(testIndex, :, :);
            % generate shuffled training data
            nSessionTrainData = shuffle3DSessionData(nSessionTrainData, trainingTargets);
            nSessionData     = [nSessionTestData; nSessionTrainData];
            shfDecodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
        end
        
        
        
        
        %% plots
        subplot(mRow, mCol, nSession)
        hold on
        
        shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(simDecodability, 1)),...
            squeeze(std(simDecodability, [], 1))/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', 'k'}, 0.5);
        
        shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(shfDecodability, 1)),...
            squeeze(std(shfDecodability, [], 1))/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', 'm'}, 0.5);
        
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([0.5 1])
        legend('boxoff')
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        xlabel('Time (s)');
        ylabel('Decodability');
        title(['# units: ' num2str(numUnits) '; # trials: ' num2str(numTrials)])
        set(gca, 'TickDir', 'out')
        
    end   
    
    setPrint(8*mCol, 6*mRow, [PlotDir 'SimultaneousUnitsDecodability/SimultaneousLDA_' DataSetList(nData).name])
end



margNames = {'Simultaneous', 'Shuffle'};
cmap = {'k', 'm'};
figure;
hold on
for nColor = 1:length(margNames)
    plot(0, nColor, 's', 'color', cmap{nColor}, 'MarkerFaceColor',cmap{nColor},'MarkerSize', 8)
    text(1, nColor, margNames{nColor})
end
xlim([0 10])
hold off
axis off
setPrint(3, 2, [PlotDir 'SimultaneousUnitsDecodability/SimultaneousLDA_Label'])

close all;

