addpath('../Func');
setDir;
load ([TempDatDir 'DataListC2SMCMCSingleTrialModel.mat']);


% for nData     = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
%     numTimeBin          = size(nDataSet(1).unit_yes_trial, 2);
%     yesProfileMatrix    = nan(length(nDataSet), numTimeBin);
%     noProfileMatrix     = yesProfileMatrix;
%     for nUnit        = 1:length(nDataSet)
%         yesData      = mean(nDataSet(nUnit).unit_yes_trial);
%         noData       = mean(nDataSet(nUnit).unit_no_trial);
%         maxData      = max([yesData, noData]);
%         minData      = min([yesData, noData]);
%         rData        = (maxData - minData);
%         yesData      = (yesData - minData)/(maxData - minData);
%         noData       = (noData - minData)/(maxData - minData);
%         yesProfileMatrix(nUnit, :)    = yesData;
%         noProfileMatrix(nUnit, :)     = noData;
%     end
%     actMat        = [yesProfileMatrix, noProfileMatrix];
%     [~, maxId]    = max(actMat, [], 2);
% 
%     timeStep  = DataSetList(nData).params.timeSeries;
%     numTime   = length(timeStep);
%     polein    = DataSetList(nData).params.polein;
%     poleout   = DataSetList(nData).params.poleout;
%     
%     
%     countMaxId = hist(maxId, 1:numTimeBin*2)/length(nDataSet)*100;
%     disp(sqrt(mean((countMaxId - 1/numTimeBin/2*100).^2)))
%     figure;
%     hold on;
%     plot(timeStep, countMaxId(1:numTimeBin), '-', 'linewid', 1.0, 'color', [0 0 0.7]);
%     plot(timeStep, countMaxId(1+numTimeBin:end), '-', 'linewid', 1.0, 'color', [0.7 0 0]);
%     xlim([timeStep(1) 2])
%     ylim([0 6])
%     gridxy ([polein, poleout, 0],[1/numTimeBin/2*100], 'Color','k','Linestyle','--','linewid', 1.0)
%     box off
%     ylabel('% Max peak')
%     xlabel('Time')
%     hold off
%     setPrint(8*2, 6, [PlotDir 'C2SModel\SingleUnitsMaxLocation_' DataSetList(nData).name '_withOLRemoval'])
% end



for nData     = [1 2]
    load([TempDatDir DataSetList(nData).name '.mat'])
    numTimeBin          = size(nDataSet(1).unit_yes_trial, 2);
    yesProfileMatrix    = nan(length(nDataSet), numTimeBin);
    noProfileMatrix     = yesProfileMatrix;
    positivePeak        = false(length(nDataSet));
    for nUnit        = 1:length(nDataSet)
        yesData      = mean(nDataSet(nUnit).unit_yes_trial);
        noData       = mean(nDataSet(nUnit).unit_no_trial);
        maxData      = max([yesData, noData]);
        minData      = min([yesData, noData]);
        rData        = (maxData - minData);
        yesData      = (yesData - minData)/(maxData - minData);
        noData       = (noData - minData)/(maxData - minData);
        yesProfileMatrix(nUnit, :)    = yesData;
        noProfileMatrix(nUnit, :)     = noData;
        positivePeak(nUnit)       = mean(yesData(1:8)) <= mean(yesData(9:47)) ...
                                   || mean(noData(1:8)) <= mean(noData(9:47));
    end
    actMat        = [yesProfileMatrix, noProfileMatrix];
    actMat        = actMat(positivePeak, :);
    [~, maxId]    = max(actMat, [], 2);

    timeStep  = DataSetList(nData).params.timeSeries;
    numTime   = length(timeStep);
    polein    = DataSetList(nData).params.polein;
    poleout   = DataSetList(nData).params.poleout;
    
    
    countMaxId = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
    disp(sqrt(mean((countMaxId - 1/numTimeBin/2*100).^2)))
    figure;
    hold on;
%     stairs(timeStep, countMaxId(1:numTimeBin), '-', 'linewid', 1.0, 'color', [0.7 0 0])
    bplot = bar(timeStep, countMaxId(1:numTimeBin), 1, 'facecolor', 'b', 'edgecolor', 'none');
    bplot.FaceAlpha = 0.5;
%     stairs(timeStep, countMaxId(1+numTimeBin:end), '-', 'linewid', 1.0, 'color', [0 0 0.7]);
    bplot = bar(timeStep, countMaxId(1+numTimeBin:end), 1, 'facecolor', 'r', 'edgecolor', 'none');
    bplot.FaceAlpha = 0.5;
    xlim([timeStep(1) 2])
    ylim([0 8])
    gridxy ([polein, poleout, 0],[1/(numTimeBin-8)/2*100], 'Color','k','Linestyle','--','linewid', 1.0)
    box off
    ylabel('% Max peak')
    xlabel('Time')
    hold off
    set(gca, 'TickDir', 'out')
    setPrint(8, 3, [PlotDir 'SingleUnitsPeakLocation\SingleUnitsMaxLocationPosNeuron_' DataSetList(nData).name])
end


close all