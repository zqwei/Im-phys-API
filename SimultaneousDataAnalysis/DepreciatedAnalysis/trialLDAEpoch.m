%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decoding which trial period one is in, how fast can you tell the data
% that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListSimultaneous.mat']);
addNoise = [1 0 0 0];
numFold  = 10;

if ~exist([PlotDir '/SimultaneousUnitsDecodabilityEpoch'],'dir')
    mkdir([PlotDir '/SimultaneousUnitsDecodabilityEpoch'])
end

for nData           = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])   
    numSession           = length(nDataSet);
    numSimultaneousUnits = zeros(numSession, 1);
    numSimultaneousTrial = zeros(numSession, 1);

    timePoints            = timePointTrialPeriod(DataSetList(nData).params.polein, DataSetList(nData).params.poleout, DataSetList(nData).params.timeSeries);
    numPeriods            = length(timePoints) - 1;
      
    
    for nSession      = 1:length(nDataSet)
        nYesDataSet   = nDataSet(nSession).unit_yes_trial;
        numYesTrial   = size(nYesDataSet, 2);
        nNoDataSet    = nDataSet(nSession).unit_no_trial;
        numNoTrial    = size(nNoDataSet, 2);
        nSessionData  = [permute(nYesDataSet, [2 1 3]); permute(nNoDataSet, [2 1 3])];
        nSessionTrial = [true(numYesTrial, 1); false(numNoTrial, 1)];
        numSimultaneousUnits(nSession) = size(nSessionData, 2);
        numSimultaneousTrial(nSession) = size(nSessionData, 1);
        decodability  = zeros(numFold, numPeriods, size(nSessionData, 3));  
        nSessionData = permute(nSessionData,[1 3 2]);
        
        for nFold        = 1:numFold
            nSessionData = nSessionData(randperm(size(nSessionData, 1)), :, :);
            EpochIndex   = epochIndex(DataSetList(nData).params);
            EpochIndex   = EpochIndex(:,ones(1,size(nSessionData, 1)))';
            numTestTrials = 10;
            decodability(nFold,:,:)    = decodabilitySliceDataTaperAccumulated(nSessionData, EpochIndex, 0, addNoise(nData), numTestTrials, numPeriods);
        end
        
        figure;
        subplot(10, 1, 1)
        [~, maxIndex] = max(mean(decodability, 1), [], 2);
        maxIndex = squeeze(maxIndex);
        imagesc(1, DataSetList(nData).params.timeSeries, maxIndex', [1 4]);
        axis off

        subplot(10, 1, 2:10)
        hold on
        area(DataSetList(nData).params.timeSeries,squeeze(mean(decodability,1))','Edgecolor','none');
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([0 1])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[0.5], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        ylabel('Decodability')
        xlabel('Time (s)')
        title(['# units: ' num2str(size(nSessionData, 3)) '; # trials: ' num2str(size(nSessionData, 1))])

        setPrint(12, 9, [PlotDir 'SimultaneousUnitsDecodabilityEpoch/SimultaneousUnitsDecodabilityEpoch_' DataSetList(nData).name '_Session_' num2str(nSession, '%02d')], 'pdf')
        
    end
    
    save([TempDatDir 'SimultaneousUnitStats_' DataSetList(nData).name '.mat'], 'numSimultaneousUnits', 'numSimultaneousTrial');
       
end

close all
