%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Different ROC
% Same number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;

load ([TempDatDir 'DataListSimultaneous.mat']);

if ~exist([PlotDir '/SimultaneousUnitsDecodability'],'dir')
    mkdir([PlotDir '/SimultaneousUnitsDecodability'])
end


% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     
%     coeffVals         = cell(length(nDataSet), 1);
%     correctRates      = cell(length(nDataSet), 1);
%     m                 = ceil(sqrt(length(nDataSet)));
%     figure;
%     
%     for nSession      = 1:length(nDataSet)
%         nYesDataSet   = nDataSet(nSession).unit_yes_trial;
%         numYesTrial   = size(nYesDataSet, 2);
%         nNoDataSet    = nDataSet(nSession).unit_no_trial;
%         numNoTrial    = size(nNoDataSet, 2);
%         nSessionData  = [permute(nYesDataSet, [2 1 3]); permute(nNoDataSet, [2 1 3])];
%         nSessionTrial = [true(numYesTrial, 1); false(numNoTrial, 1)];
%         [coeffVal, ~, ~, correctRate] = coeffSLDA(nSessionData, nSessionTrial);
%         coeffVals{nSession}           = coeffVal;
%         correctRates{nSession}        = correctRate;
%         
%         subplot(m, m, nSession);
%         plot(DataSetList(nData).params.timeSeries, correctRate);
%         xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%         ylim([0.5 1])
%         legend('boxoff')
%         gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
%         box off;
%         hold off;
%         xlabel('Time (s)');
%         ylabel('Decodability');
%         
%     end
%     
%     save([TempDatDir 'SimultaneousLDA_' DataSetList(nData).name '.mat'], 'coeffVals', 'correctRates');
%     
%     setPrint(8*m, 6*m, [PlotDir 'SimultaneousUnitsDecodability/SimultaneousUnitsDecodability_' DataSetList(nData).name], 'pdf')
% end


for nData             = 1:length(DataSetList)
    load([TempDatDir 'SimultaneousLDA_' DataSetList(nData).name '.mat'], 'correctRates');
    load([TempDatDir 'SimultaneousUnitStats_' DataSetList(nData).name '.mat'], 'numSimultaneousUnits', 'numSimultaneousTrial');
    
    [ss, sessionIndex] = sortrows([numSimultaneousUnits, numSimultaneousTrial], [-1 -2]);
    
    m                 = ceil(sqrt(length(correctRates)));
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