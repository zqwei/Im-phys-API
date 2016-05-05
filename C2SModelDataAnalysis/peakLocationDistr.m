addpath('../Func');
setDir;
load ([TempDatDir 'DataListC2SModel.mat']);


for nData     = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    numTimeBin          = size(nDataSet(1).unit_yes_trial, 2);
    yesProfileMatrix    = nan(length(nDataSet), numTimeBin);
    noProfileMatrix     = yesProfileMatrix;
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
    end
    actMat        = [yesProfileMatrix, noProfileMatrix];
    [~, maxId]    = max(actMat, [], 2);

    timeStep  = DataSetList(nData).params.timeSeries;
    numTime   = length(timeStep);
    polein    = DataSetList(nData).params.polein;
    poleout   = DataSetList(nData).params.poleout;
    
    
    countMaxId = hist(maxId, 1:numTimeBin*2)/length(nDataSet)*100;
    figure;
    hold on;
    plot(timeStep, countMaxId(1:numTimeBin), '-', 'linewid', 1.0, 'color', [0.7 0 0]);
    plot(timeStep, countMaxId(1+numTimeBin:end), '-', 'linewid', 1.0, 'color', [0 0 0.7]);
    xlim([timeStep(1) 2])
    ylim([0 6])
    gridxy ([polein, poleout, 0],[1/numTimeBin/2*100], 'Color','k','Linestyle','--','linewid', 1.0)
    box off
    ylabel('% Max peak')
    xlabel('Time')
    hold off
    setPrint(8*2, 6, [PlotDir 'C2SModel\SingleUnitsMaxLocation_' DataSetList(nData).name '_withOLRemoval'])
end

close all