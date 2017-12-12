addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

if ~exist([PlotDir 'SingleUnitsPeakLocation'],'dir')
    mkdir([PlotDir 'SingleUnitsPeakLocation'])
end


load([TempDatDir DataSetList(1).name '.mat'])
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
actMatOri     = [yesProfileMatrix, noProfileMatrix];
timeTag       = 8:60;
numFold       = 30;

sMat          = nan(numFold, 10);

for frThres = [1 4 10] % spike count in this case
    load(['validMat_' num2str(frThres, '%02d')], 'validMat')
    for nFold = 1:numFold
        actMat        = actMatOri(positivePeak & validMat(:,nFold), :);
        [~, maxId]    = max(actMat, [], 2);
        numTime       = 77;
        countMaxId    = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
        sMat(nFold, frThres)   = std(countMaxId([timeTag, timeTag+77]));
    end
end

figure;
boxplot(sMat(:, [1 4 10]));
ylim([1 2])
xlabel('Cell group')
ylabel('Peakiness')
set(gca, 'TickDir', 'out')
setPrint(8, 6, 'Peakiness', 'pdf')