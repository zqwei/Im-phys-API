%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% profileAnalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function profileAnalysisKMeans
addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CModel.mat']);
% resample all data with a larger binsize
%     timeBinLength = 3;
k             = 9;
%     m             = sqrt(k);


for nData     = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotKMeans(nDataSet, k, DataSetList(nData).params);
%         timeStep  = DataSetList(nData).params.timeSeries(ceil(timeBinLength/2):timeBinLength:end-1);
%         polein    = DataSetList(nData).params.polein;
%         poleout   = DataSetList(nData).params.poleout;
%         load([TempDatDir DataSetList(nData).name '.mat'])
%         actMat    = computeMeanActivityMatrix(nDataSet, timeBinLength);
%         idx       = kmeans(actMat, k);
%         trialSize = size(actMat, 2)/2;
%         figure;
%         for nn    = 1:k
%             subplot(m, m, nn)
%             hold on
%             if sum(idx==nn)>1
%                 mean_yes = mean(actMat(idx==nn, 1:trialSize));
%                 var_yes  = std(actMat(idx==nn, 1:trialSize))/sqrt(sum(idx==nn));
%                 shadedErrorBar(timeStep, mean_yes, var_yes, {'-b', 'linewid', 1.0}, 0.5);
%                 mean_no  = mean(actMat(idx==nn, 1+trialSize:end));
%                 var_no   = std(actMat(idx==nn, 1+trialSize:end))/sqrt(sum(idx==nn));
%                 shadedErrorBar(timeStep, mean_no, var_no, {'-r', 'linewid', 1.0}, 0.5);
%             else
%                 plot(timeStep, actMat(idx==nn, 1:trialSize), '-', 'linewid', 2.0, 'color', [0 0 0.7])
%                 plot(timeStep, actMat(idx==nn, 1+trialSize:end), '-', 'linewid', 2.0, 'color', [0.7 0 0])
%             end
%             xlabel('Normalized activity')
%             ylabel('Time')
%             xlim([timeStep(1) timeStep(end)])
%             ylim([0 1])
%             gridxy ([polein, poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
%             title(['%' num2str(sum(idx==nn)/length(idx)*100, '%.2f') ' of pop.'])
%             hold off
%             box off
%         end
    setPrint(8*ceil(sqrt(k)), 6*ceil(sqrt(k)), [PlotDir 'S2CModel/CentriodActivityProfileKMeans_K_' num2str(k, '%03d') '_' DataSetList(nData).name  '_withOLRemoval'])
end

close all
% end


% function actMat      = computeMeanActivityMatrix(nDataSet, timeBinLength)
%     numTimeBin       = floor(size(nDataSet(1).unit_yes_trial, 2)/timeBinLength);
%     trunctSize       = mod(size(nDataSet(1).unit_yes_trial, 2), timeBinLength);
%     yesProfileMatrix = nan(length(nDataSet), numTimeBin);
%     noProfileMatrix  = yesProfileMatrix;
%     
%     for nUnit        = 1:length(nDataSet)
%         yesData      = mean(nDataSet(nUnit).unit_yes_trial(:, 1:end-trunctSize));
%         noData       = mean(nDataSet(nUnit).unit_no_trial(:, 1:end-trunctSize));
%         maxData      = max([yesData, noData]);
%         minData      = min([yesData, noData]);
%         yesData      = (yesData - minData)/(maxData - minData);
%         yesData      = reshape(yesData, [timeBinLength, numTimeBin]);
%         yesProfileMatrix(nUnit, :) = mean(yesData);
% %         maxData      = max([yesData, noData]);
% %         minData      = min([yesData, noData]);
%         noData       = (noData - minData)/(maxData - minData);
%         noData       = reshape(noData, [timeBinLength, numTimeBin]);
%         noProfileMatrix(nUnit, :) = mean(noData);
%     end
%     
%     actMat  = [yesProfileMatrix, noProfileMatrix];
%     
% end

