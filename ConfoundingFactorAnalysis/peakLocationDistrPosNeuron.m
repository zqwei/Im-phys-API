addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

if ~exist([PlotDir 'SingleUnitsPeakLocation'],'dir')
    mkdir([PlotDir 'SingleUnitsPeakLocation'])
end

for nData = 10%1:length(DataSetList)
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end
    depth               = [DataSetList(nData).cellinfo(:).depth];
    depth               = depth(~neuronRemoveList)';
    
    params      = DataSetList(nData).params;
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);

    numTimeBin          = size(nDataSet(1).unit_yes_trial, 2);
    yesProfileMatrix    = nan(length(nDataSet), numTimeBin);
    noProfileMatrix     = yesProfileMatrix;
    positivePeak        = false(length(nDataSet), 1);
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
        positivePeak(nUnit)       = mean(yesData(timePoints(1):timePoints(2))) <= mean(yesData(timePoints(2):timePoints(4))) ...
                                   || mean(noData(timePoints(1):timePoints(2))) <= mean(noData(timePoints(2):timePoints(4)));

    end
    actMat        = [yesProfileMatrix, noProfileMatrix];
    depthActMat   = actMat(positivePeak & depth<200, :);
    [~, maxId]    = max(depthActMat, [], 2);
    countMaxId = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
    timeTag   = timePoints(2):timePoints(4)+13; % sample to response
    newCounts  = countMaxId([timeTag, timeTag+numTimeBin]);
    newCounts  = newCounts(newCounts>0)*numTimeBin*2/sum(newCounts>0);
    barplot(1, 1) = std(newCounts);
    [bootstat,~]  = bootstrp(1000,@std,newCounts);
    barplot(1, 2) = std(bootstat);
    disp(sum(positivePeak & depth<=200))
    bootstat_ = bootstat;
    
    depthActMat   = actMat(positivePeak & depth>200 & depth<400, :);
    [~, maxId]    = max(depthActMat, [], 2);
    countMaxId = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
    newCounts  = countMaxId([timeTag, timeTag+numTimeBin]);
    newCounts  = newCounts(newCounts>0)*numTimeBin*2/sum(newCounts>0);
    barplot(2, 1) = std(newCounts);
    [bootstat,~]  = bootstrp(1000,@std,newCounts);
    barplot(2, 2) = std(bootstat);
    disp(sum(positivePeak & depth>200 & depth<400))
    
    [h, p] = ttest2(bootstat_, bootstat);
    disp(p)
    
    figure;
    hold on;
    bar(1:2, barplot(:, 1));
    errorbar(1:2, barplot(:,1), barplot(:, 2));
    ylabel('Peakiness')
    xlabel('Depth (um)')
    ylim([0 1.6])
    xlim([0.5 2.5])
    set(gca, 'TickDir', 'out')
    set(gca, 'XTick', 1:2, 'YTick', [0 1.6])
    setPrint(8, 6, ['Peakiness_' DataSetList(nData).name], 'pdf')
end

close all