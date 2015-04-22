%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_8

addpath('../Func');
setDir;
load ([TempDatDir 'DataList.mat']);
addNoise         = [1 0 0 0 0 0];



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 8.1  Decoding which trial period one is in, how fast can you tell the data
% %     that the trial period switched
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])   
%     m                 = ceil(sqrt(length(nDataSet3D)));
%     figure;
%     for nSession = 1:length(nDataSet3D)
%         nSessionData = [permute(nDataSet3D(nSession).unit_yes_trial,[2 1 3]); permute(nDataSet3D(nSession).unit_no_trial,[2 1 3])];
%         nSessionData = permute(nSessionData,[1 3 2]);
%         numTrials    = size(nSessionData ,1);
%         EpochIndex   = epochIndex(DataSetList(nData).params);
%         EpochIndex   = EpochIndex(:,ones(1,numTrials))';
%         numTestTrials = round(numTrials*0.2);
%         decodability = decodabilitySliceDataTaper(nSessionData, EpochIndex, 0, addNoise(nData), numTestTrials);
%         subplot(m, m, nSession);
%         hold on
%         plot(DataSetList(nData).params.timeSeries, decodability, 'k', 'linewid',1);
%         xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%         ylim([0 1])
%         gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
%         box off;
%         hold off;
%         ylabel('Decodability')
%         xlabel('Time (s)')
%     end
%     setPrint(6*m, 4.5*m, [PlotDir '/Collected_Units_Decodability_EpochLDA1_' DataSetList(nData).name], 'pdf')
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 8.2  Decoding which trial period one is in, how fast can you tell the data
% %     that the trial period switched (Multi-taper)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])   
%     m                 = ceil(sqrt(length(nDataSet3D)));
%     figure;
%     for nSession = 1:length(nDataSet3D)
%         nSessionData = [permute(nDataSet3D(nSession).unit_yes_trial,[2 1 3]); permute(nDataSet3D(nSession).unit_no_trial,[2 1 3])];
%         nSessionData = permute(nSessionData,[1 3 2]);
%         numTrials    = size(nSessionData ,1);
%         EpochIndex   = epochIndex(DataSetList(nData).params);
%         EpochIndex   = EpochIndex(:,ones(1,numTrials))';
%         numTestTrials = round(numTrials*0.2);
%         numTaper     = 9;
%         subplot(m, m, nSession);
%         hold on
%         for nTaper   = 0:numTaper
%             nColor       = [nTaper, 0, 0] * 0.1;
%             decodability = decodabilitySliceDataTaper(nSessionData, EpochIndex, nTaper, addNoise(nData), numTestTrials);
%             plot(DataSetList(nData).params.timeSeries(1:end-numTaper), decodability(1:end+nTaper-numTaper), 'Color', nColor, 'linewid',1);
%         end  
%         xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%         ylim([0 1])
%         gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
%         box off;
%         hold off;
%         ylabel('Decodability')
%         xlabel('Time (s)')
%     end
%     setPrint(6*m, 4.5*m, [PlotDir '/Collected_Units_Decodability_EpochLDATaper1_' DataSetList(nData).name], 'pdf')
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.3  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched (Multi-LDA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])   
%     m                 = ceil(sqrt(length(nDataSet3D)));
%     numFold           = 10;
%     figure;
%     for nSession     = 1:length(nDataSet3D)
%         nSessionData = [permute(nDataSet3D(nSession).unit_yes_trial,[2 1 3]); permute(nDataSet3D(nSession).unit_no_trial,[2 1 3])];
%         % nSessionData : Ntrial x Nneuron x Nt
%         nSessionData = permute(nSessionData,[1 3 2]);
%         % nSessionData : Ntrial x Nt x Nneuron
%         numTrials    = size(nSessionData ,1);
%         EpochIndex   = epochIndex(DataSetList(nData).params);
%         EpochIndex   = EpochIndex(:,ones(1,numTrials))';
%         % EpochIndex : Ntrial x Nt
%         
%         numTestTrials = round(numTrials*0.2);
%         subplot(m, m, nSession);
%         hold on
%         timePoints   = timePointTrialPeriod(DataSetList(nData).params.polein, DataSetList(nData).params.poleout, DataSetList(nData).params.timeSeries);
%         for nPeriods        = 1: length(timePoints) -2 
%             nColor          = [nPeriods, 0, 0] * 0.3;
%             sliceIndex      = timePoints(nPeriods)+3:timePoints(nPeriods+2)-3;% timePoints(nPeriods+1)-7:timePoints(nPeriods+1)+6;% timePoints(nPeriods)+3:timePoints(nPeriods+2)-3;
%             decodability    = zeros(numFold, length(sliceIndex));
%             for nFold       = 1:numFold
%                 decodability(nFold,:)    = decodabilitySliceDataPeriod(nSessionData, EpochIndex, sliceIndex, addNoise(nData), numTestTrials);
%             end
%             % plot(DataSetList(nData).params.timeSeries(sliceIndex), mean(decodability), 'Color', nColor, 'linewid',1);
%             shadedErrorBar(DataSetList(nData).params.timeSeries(sliceIndex), mean(decodability), std(decodability), {'Color', nColor, 'linewid',1});
%         end        
%         xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%         ylim([0 1])
%         gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[0.5], 'Color','k','Linestyle','--','linewid', 0.5)
%         box off;
%         hold off;
%         ylabel('Decodability')
%         xlabel('Time (s)')
%     end
%     setPrint(6*m, 4.5*m, [PlotDir '/Collected_Units_Decodability_EpochLDA3_' DataSetList(nData).name], 'pdf')
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 8.4  Decoding which trial period one is in, how fast can you tell the data
% %     that the trial period switched (Multi-LDA) Based on 8.3 New plot
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])   
%     m                 = ceil(length(nDataSet3D)/4);
%     numFold           = 10;
%     figure;
%     
%     numUnit          = arrayfun(@(x) length(x.nUnit), nDataSet3D, 'UniformOutput', false);
%     [~, sortUnit]    = sort([numUnit{:}],'descend');
%     
%     nDataSet3D       = nDataSet3D(sortUnit);
%     
%     for nSession     = 1:length(nDataSet3D)
%         nSessionData = [permute(nDataSet3D(nSession).unit_yes_trial,[2 1 3]); permute(nDataSet3D(nSession).unit_no_trial,[2 1 3])];
%         % nSessionData : Ntrial x Nneuron x Nt
%         nSessionData = permute(nSessionData,[1 3 2]);
%         % nSessionData : Ntrial x Nt x Nneuron
%         numTrials    = size(nSessionData ,1);
%         EpochIndex   = epochIndex(DataSetList(nData).params);
%         EpochIndex   = EpochIndex(:,ones(1,numTrials))';
%         % EpochIndex : Ntrial x Nt
%         
%         numTestTrials = round(numTrials*0.2);
%         subplot(m, 4, nSession);
%         hold on
%         timePoints   = timePointTrialPeriod(DataSetList(nData).params.polein, DataSetList(nData).params.poleout, DataSetList(nData).params.timeSeries);
%         numPeriods   = length(timePoints) - 1;
%         for nPeriods        = 1: length(timePoints) -2 
%             nColor          = [nPeriods, 0, 0] * 0.3;
%             sliceIndex      = timePoints(nPeriods)+3:timePoints(nPeriods+2)-3;% timePoints(nPeriods+1)-7:timePoints(nPeriods+1)+6;% timePoints(nPeriods)+3:timePoints(nPeriods+2)-3;
%             decodability    = zeros(numFold, numPeriods, length(sliceIndex));
%             for nFold       = 1:numFold
%                 decodability(nFold,:,:)    = decodabilitySliceDataPeriodAccumulated(nSessionData, EpochIndex, sliceIndex, addNoise(nData), numTestTrials, numPeriods);
%             end
%             area(DataSetList(nData).params.timeSeries(sliceIndex),squeeze(mean(decodability,1))','Edgecolor','none');
%         end        
%         xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%         ylim([0 1])
%         gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[0.5], 'Color','k','Linestyle','--','linewid', 0.5)
%         box off;
%         hold off;
%         ylabel('Decodability')
%         xlabel('Time (s)')
%         title([num2str(size(nSessionData, 3)) ' Neurons'])
%         ylim([0 1])
%         xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     end
%     setPrint(4*4, 3*m, [PlotDir '/Collected_Units_Decodability_EpochLDA3Acc_' DataSetList(nData).name], 'pdf')
% end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 8.5  Decoding which trial period one is in, how fast can you tell the data
% %     that the trial period switched (Single-LDA) Based on 8.1 New plot
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numFold               = 10;
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     numUnit          = arrayfun(@(x) length(x.nUnit), nDataSet3D, 'UniformOutput', false);
%     [~, sortUnit]    = sort([numUnit{:}],'descend');
%     nDataSet3D       = nDataSet3D(sortUnit);
%     m                 = ceil(length(nDataSet3D)/4);
%     figure;
%     for nSession = 1:length(nDataSet3D)
%         nSessionData = [permute(nDataSet3D(nSession).unit_yes_trial,[2 1 3]); permute(nDataSet3D(nSession).unit_no_trial,[2 1 3])];
%         % nSessionData : Ntrial x Nneuron x Nt
%         nSessionData = permute(nSessionData,[1 3 2]);
%         % nSessionData : Ntrial x Nt x Nneuron
%         numTrials    = size(nSessionData ,1);
%         EpochIndex   = epochIndex(DataSetList(nData).params);
%         EpochIndex   = EpochIndex(:,ones(1,numTrials))';
%         % EpochIndex : Ntrial x Nt
%         numTestTrials = round(numTrials*0.2);
%         timePoints   = timePointTrialPeriod(DataSetList(nData).params.polein, DataSetList(nData).params.poleout, DataSetList(nData).params.timeSeries);
%         numPeriods   = length(timePoints) - 1;
%         decodability = zeros(numFold, numPeriods, size(nSessionData,2));
%         for nFold       = 1:numFold
%             decodability(nFold,:,:)    = decodabilitySliceDataTaperAccumulated(nSessionData, EpochIndex, 0, addNoise(nData), numTestTrials, numPeriods);
%         end
%         subplot(m, 4, nSession);
%         hold on
%         area(DataSetList(nData).params.timeSeries,squeeze(mean(decodability,1))','Edgecolor','none');
%         xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%         ylim([0 1])
%         gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[0.5], 'Color','k','Linestyle','--','linewid', 0.5)
%         box off;
%         hold off;
%         ylabel('Decodability')
%         xlabel('Time (s)')
%         title([num2str(size(nSessionData, 3)) ' Neurons'])
%     end
%     setPrint(4*4, 3*m, [PlotDir '/Collected_Units_Decodability_EpochLDA1ACC_' DataSetList(nData).name], 'pdf')
% end

if ~exist([PlotDir '/Collected_Units_Decodability_Epoch'],'dir')
    mkdir([PlotDir '/Collected_Units_Decodability_Epoch'])
end
controlUnits        = 35;
numTrials           = 400;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
trainingTargets     = trainingTargets(randperm(numTrainingTrials));
testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
testTargets         = testTargets(randperm(numTestTrials));
totTargets          = [testTargets; trainingTargets];

numFold               = 30;

for nData             = [1 3 4 5 6]%1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])   
    figure;
    timePoints   = timePointTrialPeriod(DataSetList(nData).params.polein, DataSetList(nData).params.poleout, DataSetList(nData).params.timeSeries);
    numPeriods   = length(timePoints) - 1;
    decodability = zeros(numFold, numPeriods, size(nDataSet(1).unit_yes_trial,2));
    for nFold       = 1:numFold
        numUnits     = length(nDataSet);
        nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, controlUnits)), totTargets, numTestTrials);
        % nSessionData : Ntrial x Nneuron x Nt
        nSessionData = permute(nSessionData,[1 3 2]);
        % nSessionData : Ntrial x Nt x Nneuron
        % numTrials    = size(nSessionData ,1);
        EpochIndex   = epochIndex(DataSetList(nData).params);
        EpochIndex   = EpochIndex(:,ones(1,numTrials))';
        % EpochIndex : Ntrial x Nt
        % numTestTrials = round(numTrials*0.2);        
        decodability(nFold,:,:)    = decodabilitySliceDataTaperAccumulated(nSessionData, EpochIndex, 0, addNoise(nData), numTestTrials, numPeriods);
    end
    hold on
    area(DataSetList(nData).params.timeSeries,squeeze(mean(decodability,1))','Edgecolor','none');
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0 1])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[0.5], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    ylabel('Decodability')
    xlabel('Time (s)')
    setPrint(4, 3, [PlotDir 'Collected_Units_Decodability_Epoch/Collected_Units_Decodability_EpochLDA1AllAccFixedNumUnits_' DataSetList(nData).name], 'pdf')
end

close all