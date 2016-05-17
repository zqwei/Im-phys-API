%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Same ROC
% Different number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('../Func');
setDir;


load ([TempDatDir 'DataListShuffle.mat']);

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


stepSize            = 10;
numFold             = 10;
thres               = 0.65;

figure;
hold on
for nData             = [1 3 4]
    load(['decodabilityAll_' DataSetList(nData).name], 'decodabilityAll')
    maxRandPickUnits  = 20;
    delayMean         = zeros(maxRandPickUnits, 1);
    delayStd          = zeros(maxRandPickUnits, 1);
    
    for numRandPickUnits      = 1:maxRandPickUnits;
        delays        = zeros(numFold, 1);
        for nFold     = 1:numFold
            delayIndex = find(squeeze(decodabilityAll(numRandPickUnits, nFold, :))>thres);
            delayIndex = delayIndex(delayIndex>8);
            delays(nFold) = DataSetList(nData).params.timeSeries(min(delayIndex)); 
        end
        delayMean(numRandPickUnits) = mean(delays)-(-2.6);
        delayStd(numRandPickUnits)  = sem(delays);
    end
    delayMean = delayMean * 1000;
    delayStd  = delayStd * 1000;
    
    
    shadedErrorBar((1:maxRandPickUnits)*stepSize, delayMean, delayStd,...
        {'-', 'linewid', 1.0, 'color', cmap(nData,:)}, 0.5);  
end

xlim([0 maxRandPickUnits*stepSize+5]);
% ylim([0.5 1])
set(gca, 'TickDir', 'out')
box off;
hold off;
xlabel('# Units in analysis');
ylabel('Time > 0.65 (ms)');

setPrint(8, 6, [PlotDir 'CollectedUnitsDecodability/CollectedUnitsDecodabilityROC_DelayInSample'])


figure;
hold on
for nData             = [1 3 4]
    load(['decodabilityAll_' DataSetList(nData).name], 'decodabilityAll')
    maxRandPickUnits  = 20;
    accuracyMean         = zeros(maxRandPickUnits, 1);
    accuracyStd          = zeros(maxRandPickUnits, 1);
    
    for numRandPickUnits      = 1:maxRandPickUnits;
        accuracys        = zeros(numFold, 1);
        for nFold     = 1:numFold
            accuracys(nFold) = mean(squeeze(decodabilityAll(numRandPickUnits, nFold, 47:end))); 
        end
        accuracyMean(numRandPickUnits) = mean(accuracys);
        accuracyStd(numRandPickUnits)  = sem(accuracys);
    end
    
    shadedErrorBar((1:maxRandPickUnits)*stepSize, accuracyMean, accuracyStd,...
        {'-', 'linewid', 1.0, 'color', cmap(nData,:)}, 0.5);  
end

xlim([0 maxRandPickUnits*stepSize+5]);
ylim([0.7 1])
set(gca, 'TickDir', 'out')
box off;
hold off;
xlabel('# Units in analysis');
ylabel('Mean accuracy in response epoch');
setPrint(8, 6, [PlotDir 'CollectedUnitsDecodability/CollectedUnitsDecodabilityROC_MeanPerformanceInResponse'])
