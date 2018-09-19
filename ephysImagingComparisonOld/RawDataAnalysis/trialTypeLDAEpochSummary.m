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

stepSize            = 10;
numFold             = 10;


for nData             = [1 3]
    load(['decodabilityAll_' DataSetList(nData).name], 'decodabilityAll')
    
    errorset          = nan(1, 3);
    barset            = 1:3;
    barmean           = nan(1, 3);
    
    params            = DataSetList(nData).params;
    
    timePoints        = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    
    for numRandPickUnits      = [1 5]
        for nTime     = 1:3
            barmean(nTime)    = mean(mean(decodabilityAll(numRandPickUnits, :, timePoints(nTime+1):timePoints(nTime+2)), 3));
            errorset(nTime)   = std(mean(decodabilityAll(numRandPickUnits, :, timePoints(nTime+1):timePoints(nTime+2)), 3))/sqrt(numFold);
        end
        
    
        figure;
        hold on
        errorbar(barset, barmean, errorset, '.k')
        bar(barset, barmean)
        xlim([0.5 3.5]);
        ylim([0.5 1.1]);
        set(gca, 'TickDir', 'out')
        box off;
        hold off;
        setPrint(8, 6, [PlotDir 'CollectedUnitsDecodability/CollectedUnitsDecodabilityEpochBar_idx_' num2str(numRandPickUnits) '_' DataSetList(nData).name])
    end
end