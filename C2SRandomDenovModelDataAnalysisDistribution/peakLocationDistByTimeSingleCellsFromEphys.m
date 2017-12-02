%
% Compute distribution of selectivity for each single neuron
% 
% Here we drop the test of parameter of tau_r
% 
% 
% -------------------------------------------------------------------------
% version 1.0
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

% all data is precomputed by code meanSynthEphysTraces.m

indexDatasets = [3, 4];
% indexDatasets = [2, 3, 4, 5, 6];
% 2: short Ca GP517
% 3: short Ca slow
% 4: short Ca slow virus
% 5: long Ca fast
% 6: long Ca slow

params           = DataSetList(1).params;
load([TempDatDir DataSetList(1).name '.mat'])
spikeDataSet     = nDataSet;
timeTag          = 8:60;

numTimeBin       = size(spikeDataSet(1).unit_yes_trial, 2);
yesProfileMatrix = nan(length(spikeDataSet), numTimeBin);
noProfileMatrix  = yesProfileMatrix;
positivePeak     = false(length(spikeDataSet));
for nUnit        = 1:length(spikeDataSet)
    yesUnitData      = spikeDataSet(nUnit).unit_yes_trial;
    noUnitData       = spikeDataSet(nUnit).unit_no_trial;
    yesData          = mean(yesUnitData);
    noData           = mean(noUnitData);
    maxData          = max([yesData, noData]);
    minData          = min([yesData, noData]);
    rData            = (maxData - minData);
    yesData          = (yesData - minData)/(maxData - minData);
    noData           = (noData - minData)/(maxData - minData);
    yesProfileMatrix(nUnit, :)    = yesData;
    noProfileMatrix(nUnit, :)     = noData;
    positivePeak(nUnit)       = mean(yesData(1:8)) <= mean(yesData(9:47)) ...
                               || mean(noData(1:8)) <= mean(noData(9:47));
end

actMat        = [yesProfileMatrix, noProfileMatrix];
actMat        = actMat(positivePeak, :);
[~, maxId]    = max(actMat, [], 2);
% % countMaxId    = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
% % peakiness(nFF) = std(countMaxId([timeTag, timeTag+77]));

for nData     = indexDatasets   
    load([TempDatDir 'directDeconv/' DataSetList(nData).name '_ephysSimSumMat.mat'], 'ephysSimMatSum');
    size(ephysSimMatSum)
    timeWeightVec  = zeros(numTimeBin * 2, 1);
    figure
    plot(ephysWeightVec)
    ephysWeightVec = mean(ephysSimMatSum(positivePeak, :)/97/97, 2);
    for nUnit      = 1:length(maxId)
        timeWeightVec(maxId(nUnit)) = timeWeightVec(maxId(nUnit)) + ephysWeightVec(nUnit);
    end
    plot(timeWeightVec)
    countMaxId    = timeWeightVec/sum(timeWeightVec)*100;
    peakiness     = std(countMaxId([timeTag, timeTag+77]))
%     figure
%     bar(timeWeightVec)
end

