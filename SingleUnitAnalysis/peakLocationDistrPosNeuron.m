addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

if ~exist([PlotDir 'SingleUnitsPeakLocation'],'dir')
    mkdir([PlotDir 'SingleUnitsPeakLocation'])
end

for nData     = 1:length(DataSetList)
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end
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
    timeTag   = 8:60; % sample to response
    numTime   = length(timeTag);
    polein    = DataSetList(nData).params.polein;
    poleout   = DataSetList(nData).params.poleout;
    
    countMaxId = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
    % mean(countMaxId([8:77, 95:end]))
    % disp(sqrt(mean((countMaxId([8:77, 95:end]) - 1/numTimeBin/2*100).^2)))
    std(countMaxId([timeTag, timeTag+77]))%/mean(countMaxId([timeTag, timeTag+77]))
    [bootstat,bootsam] = bootstrp(1000,@std,countMaxId([timeTag, timeTag+77]));
%     mean(bootstat)
    std(bootstat)
%     [h, p] = ttest(bootstat)
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
    set(gca,'TickDir','out')
    ylabel('% Max peak')
    xlabel('Time')
    hold off
    setPrint(8, 3, [PlotDir 'SingleUnitsPeakLocation\SingleUnitsMaxLocationPosNeuron_' DataSetList(nData).name])
    setPrint(8, 3, [PlotDir 'SingleUnitsPeakLocation\' DataSetList(nData).name '_peakness'], 'svg')
end

close all