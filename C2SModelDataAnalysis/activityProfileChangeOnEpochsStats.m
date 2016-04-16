%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change of neuronal activity as a function of behavioral epoch onset time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function activityProfileChangeOnEpochsStats
    addpath('../Func');
    setDir;
    load ([TempDatDir 'DataListC2SModel.mat']);
    smoothedBinsize = 11;
    stdThres      = 3;
    ttestThres    = 0.1;
    spThres       = 0.2;
    for nData     = 1:length(DataSetList)
        load([TempDatDir DataSetList(nData).name '.mat'])
        [actMat, splineActMat, stdActMat] = smoothedMeanActivityMatrix(nDataSet, smoothedBinsize);
        peakMat   = false(size(actMat));
        timeStep  = DataSetList(nData).params.timeSeries;
        polein    = DataSetList(nData).params.polein;
        poleout   = DataSetList(nData).params.poleout;
        numTime   = length(timeStep);
        for nUnit = 1:length(nDataSet)
            validIndex = (splineActMat(nUnit,:)>actMat(nUnit,:)+stdActMat(nUnit,:)*stdThres) | (splineActMat(nUnit,:)<actMat(nUnit,:)-stdActMat(nUnit,:)*stdThres);
            validIndex = validIndex & ~(actMat(nUnit,:)<0.2 & splineActMat(nUnit,:)<0.3) & [timeStep>polein, timeStep>polein];
            if sum(validIndex) && std(splineActMat(nUnit,:))>spThres && ttest(stdActMat(nUnit,:), ttestThres, 'tail', 'left')
                peakMat(nUnit, :) = validIndex;
            end
        end
        
        figure
        hold on
        plot(timeStep, mean(peakMat(:, 1:numTime))*100, '-r', 'linewid', 1.0)
        plot(timeStep, mean(peakMat(:, numTime+1:end))*100, '-b', 'linewid', 1.0)
        ylabel('% neurons')
        xlabel('Time')
        ylim([0 2.5])
        xlim([timeStep(1) timeStep(end)])
        gridxy ([polein, poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
        box off
        setPrint(8*2, 6, [PlotDir 'C2SModel/SingleUnitsPeakFinderStats_' DataSetList(nData).name])
    end
end

function [actMat, splineActMat, stdActMat]= smoothedMeanActivityMatrix(nDataSet, smoothedBinsize)
    numTimeBin          = size(nDataSet(1).unit_yes_trial, 2);
    yesProfileMatrix    = nan(length(nDataSet), numTimeBin);
    noProfileMatrix     = yesProfileMatrix;
    spYesProfileMatrix  = nan(length(nDataSet), numTimeBin);
    spNoProfileMatrix   = yesProfileMatrix;
    stdYesProfileMatrix = nan(length(nDataSet), numTimeBin);
    stdNoProfileMatrix  = yesProfileMatrix;
    
    for nUnit        = 1:length(nDataSet)
        yesData      = mean(nDataSet(nUnit).unit_yes_trial);
        noData       = mean(nDataSet(nUnit).unit_no_trial);
        maxData      = max([yesData, noData]);
        minData      = min([yesData, noData]);
        rData        = (maxData - minData);
        yesData      = (yesData - minData)/(maxData - minData);
        noData       = (noData - minData)/(maxData - minData);
        stdYesData   = std(nDataSet(nUnit).unit_yes_trial)/rData/sqrt(size(nDataSet(nUnit).unit_yes_trial, 1));
        stdNoData    = std(nDataSet(nUnit).unit_no_trial)/rData/sqrt(size(nDataSet(nUnit).unit_no_trial, 1));
        stdYesProfileMatrix(nUnit, :) = stdYesData;
        stdNoProfileMatrix(nUnit, :)  = stdNoData;
        yesProfileMatrix(nUnit, :)    = yesData;
        noProfileMatrix(nUnit, :)     = noData;
        yesData      = smooth(yesData, smoothedBinsize, 'moving');
        noData       = smooth(noData, smoothedBinsize, 'moving');
        spYesProfileMatrix(nUnit, :)  = yesData;
        spNoProfileMatrix(nUnit, :)   = noData;
    end
    
    actMat        = [yesProfileMatrix, noProfileMatrix];
    splineActMat  = [spYesProfileMatrix, spNoProfileMatrix];
    stdActMat     = [stdYesProfileMatrix, stdNoProfileMatrix];
end

