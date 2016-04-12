%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change of neuronal activity as a function of behavioral epoch onset time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function activityProfileChangeOnEpochs
    addpath('../Func');
    setDir;
    load ([TempDatDir 'DataListShuffle.mat']);
    smoothedBinsize = 11;
    for nData     = [1 3 4]
        load([TempDatDir DataSetList(nData).name '.mat'])
        [actMat, splineActMat, stdActMat] = smoothedMeanActivityMatrix(nDataSet, smoothedBinsize);
%         residual_mat = (actMat - splineActMat)./splineActMat;
        timeStep  = DataSetList(nData).params.timeSeries;
        numTime   = length(timeStep);
        polein    = DataSetList(nData).params.polein;
        poleout   = DataSetList(nData).params.poleout;
        for nUnit = 1:length(nDataSet)
            validIndex = (splineActMat(nUnit,:)>actMat(nUnit,:)+stdActMat(nUnit,:)*3) | (splineActMat(nUnit,:)<actMat(nUnit,:)-stdActMat(nUnit,:)*3);
            validIndex = validIndex & ~(actMat(nUnit,:)<0.2 & splineActMat(nUnit,:)<0.3) & [timeStep>polein, timeStep>polein];
            if sum(validIndex) && std(splineActMat(nUnit,:))>0.2 && ttest(stdActMat(nUnit,:), 0.1, 'tail', 'left')
                disp(std(splineActMat(nUnit,:)))
                figure
                subplot(1, 2, 1)
                hold on
                shadedErrorBar(timeStep, actMat(nUnit, 1:numTime), stdActMat(nUnit, 1:numTime)*1.5, {'-r', 'linewid', 1.0}, 0.5)
                plot(timeStep, splineActMat(nUnit, 1:numTime), '--r', 'linewid', 1.0)
                locs = validIndex(1:numTime); %abs(residual_mat(nUnit, 1:numTime))>1.0 & 
                locs = find(locs);
                if ~isempty(locs)
                    plot(timeStep(locs), actMat(nUnit, locs), 'vr', 'linewid', 1.0, 'markerfacecolor', 'r')
                end
                ylim([0 1.05])
                xlim([timeStep(1) timeStep(end)])
                ylabel('Normalized activity')
                xlabel('Time')
                gridxy ([polein, poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
                box off
                subplot(1, 2, 2)
                hold on
                shadedErrorBar(timeStep, actMat(nUnit, numTime+1:end), stdActMat(nUnit, numTime+1:end)*1.5, {'-b', 'linewid', 1.0}, 0.5)
                plot(timeStep, splineActMat(nUnit, numTime+1:end), '--b', 'linewid', 1.0)
                locs = validIndex(numTime+1:end); % abs(residual_mat(nUnit, numTime+1:end))>1.0 & 
                locs = find(locs);
                if ~isempty(locs)
                    plot(timeStep(locs), actMat(nUnit, numTime+locs), 'vb', 'linewid', 1.0, 'markerfacecolor', 'b')
                end
                ylabel('Normalized activity')
                xlabel('Time')
                ylim([0 1.05])
                xlim([timeStep(1) timeStep(end)])
                gridxy ([polein, poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
                box off
                setPrint(8*2, 6, [PlotDir 'SingleUnitsPeakFinder\' DataSetList(nData).name '_NeuronIdx_' num2str(nUnit, '%04d')])
                close all
            end
        end
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

