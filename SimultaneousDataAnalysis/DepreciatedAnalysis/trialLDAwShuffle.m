%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Different ROC
% Same number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;

load ([TempDatDir 'DataListSimultaneous.mat']);

if ~exist([PlotDir '/SimultaneousShuffleUnitsDecodability'],'dir')
    mkdir([PlotDir '/SimultaneousShuffleUnitsDecodability'])
end

numFold              = 10;

% 
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     load([TempDatDir 'SimultaneousUnitStats_' DataSetList(nData).name '.mat'], 'numSimultaneousUnits', 'numSimultaneousTrial');
%     
%     ssToAnalysis       = numSimultaneousUnits>10 & (numSimultaneousTrial> 2* numSimultaneousUnits);
%     nDataSet           = nDataSet(ssToAnalysis);
%     ss                 = [numSimultaneousUnits(ssToAnalysis), numSimultaneousTrial(ssToAnalysis)];
%     
%     coeffVals         = cell(length(nDataSet), 1+numFold);
%     correctRates      = cell(length(nDataSet), 1+numFold);
%     
%     for nSession      = 1:length(nDataSet)
%         nYesDataSet   = nDataSet(nSession).unit_yes_trial;
%         numYesTrial   = size(nYesDataSet, 2);
%         nNoDataSet    = nDataSet(nSession).unit_no_trial;
%         numNoTrial    = size(nNoDataSet, 2);
%         nSessionData  = [permute(nYesDataSet, [2 1 3]); permute(nNoDataSet, [2 1 3])];
%         numNeuron     = size(nSessionData, 2);
%         nSessionTrial = [true(numYesTrial, 1); false(numNoTrial, 1)];
%         [coeffVal, ~, ~, correctRate] = coeffSLDA(nSessionData, nSessionTrial);
%         coeffVals{nSession, 1}           = coeffVal;
%         correctRates{nSession, 1}        = correctRate;
%         tSessionData   = nSessionData;
%         
%         for nFold     = 1:numFold
%             randYesIndex   = randi(numYesTrial, numYesTrial, numNeuron);
%             randNoIndex    = randi(numNoTrial, numNoTrial, numNeuron) + numYesTrial;
%             randIndex      = [randYesIndex; randNoIndex];
%             for nUnit      = 1:numNeuron
%                 tSessionData(:, nUnit, :) = nSessionData(randIndex(:, nUnit), nUnit, :);
%             end
%             [coeffVal, ~, ~, correctRate] = coeffSLDA(tSessionData, nSessionTrial);
%             coeffVals{nSession, 1+nFold}           = coeffVal;
%             correctRates{nSession, 1+nFold}        = correctRate;
%         end
%         
%     end
%     
%     save([TempDatDir 'SimultaneousShuffleLDA_' DataSetList(nData).name '.mat'], 'coeffVals', 'correctRates', 'ss');
%     
% end


% for nData             = 1:length(DataSetList)
%     load([TempDatDir 'SimultaneousShuffleLDA_' DataSetList(nData).name '.mat'], 'correctRates', 'ss');
%     
%     for nSession      = 1:size(correctRates, 1)
%         correctRate   = cell2mat(correctRates(nSession,:)');   
%         figure;
%         hold on;
%         plot(DataSetList(nData).params.timeSeries, correctRate(2:numFold+1, :), 'Color',[0.5 0.5 0.5],'Linestyle','-','linewid', 0.5);
%         plot(DataSetList(nData).params.timeSeries, mean(correctRate(2:numFold+1, :)), 'Color','k','Linestyle','--','linewid', 1.0);
%         plot(DataSetList(nData).params.timeSeries, correctRate(1, :), 'Color','k','Linestyle','-','linewid', 1.0);
%         xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%         ylim([0.5 1])
%         legend('boxoff')
%         gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
%         box off;
%         hold off;
%         xlabel('Time (s)');
%         ylabel('Decodability');
%         setPrint(8, 6, [PlotDir 'SimultaneousShuffleUnitsDecodability/SimultaneousUnitsDecodability_' DataSetList(nData).name '_Session_' num2str(nSession)], 'pdf')
%     end
%     
% end



markerSet             = {'o', 's', 'v', '^'};

figure;
hold on
for nData             = 1:length(DataSetList)
    load([TempDatDir 'SimultaneousShuffleLDA_' DataSetList(nData).name '.mat'], 'correctRates', 'ss');
    trialVar          = zeros(size(correctRates, 1), 1);
    trialEV           = zeros(size(correctRates, 1), 1);
    for nSession      = 1:size(correctRates, 1)
        correctRate   = cell2mat(correctRates(nSession,:)');   
        trialEV(nSession)  = sqrt(mean((correctRate(1, :) - mean(correctRate(2:numFold+1, :))).^2));
        trialVar(nSession) = sqrt(mean(var(correctRate(2:numFold+1, :))));
    end
    scatter(trialVar, trialEV, ss(:,1), ss(:,2), markerSet{nData}, 'filled')
end
hold off
colormap(jet)
xlabel('Trial by trial var.')
ylabel('Estimation error')
setPrint(8, 6, [PlotDir 'SimultaneousShuffleUnitsDecodability/SimultaneousUnitsDecodability_' DataSetList(nData).name '_AllSession'], 'pdf')

close all;