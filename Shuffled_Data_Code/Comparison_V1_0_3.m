% 
% Comparison based on single unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0.3
% 
% 
% 
% Check all neurons stastical properity comparing to ephys
% 
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);


fileList                 = {SpikeFileList; CaImagingShortDelayFastFileList; CaImagingShortDelaySlowFileList;...
                            CaImagingLongDelayFastFileList; CaImagingLongDelaySlowFileList; CaImagingShortDelaySlowVirusFileList};

nData                    = 1;
disp(DataSetList(nData).name)
load([TempDatDir DataSetList(nData).name '.mat']);
ephysCellIndex           = [DataSetList(nData).cellinfo(:).depth]  > 120 &...
                           [DataSetList(nData).cellinfo(:).depth]  < 740;
% ephysCellDepth  = [DataSetList(nData).cellinfo(ephysCellIndex).depth];
% ephysCellML     = [DataSetList(nData).cellinfo(ephysCellIndex).ML_axis];
% [E, ephysCellIndex] = sortrows([ephysCellDepth', ephysCellML'], [2, 1]);


numUnits                  = length(nDataSet(ephysCellIndex));
zScores                   = applyFuncToCompareTrialType(nDataSet(ephysCellIndex), @dValue);
zScores(isnan(zScores))   = 0;     

actMat                    = abs(zScores);
numT                      = size(zScores,2);

bumpActThres              = 0.4; % > bumpActThres considering as a bump
bumpMat                   = actMat > bumpActThres;
bumpSize                  = ones(size(actMat,1),1);
bumpStartPoint            = nan(size(actMat,1),1);
bumpSign                  = ones(size(actMat,1),1);

for nUnit = 1: size(actMat,1)
    diffActMat            = diff([bumpMat(nUnit,:),0]);
    beginPoint            = find(diffActMat==1);
    endPoint              = find(diffActMat==-1);
    if length(endPoint)>length(beginPoint); endPoint = endPoint(2:end); end
    [bumplength, bumpIndex] = max(endPoint - beginPoint);
    if isempty(bumplength); bumplength = 0; end
    bumpSize(nUnit)       = bumplength;
    if bumplength==0
        bumpStartPoint(nUnit) = numT;
        bumpSign(nUnit)       = -1;
    else 
        bumpStartPoint(nUnit) = beginPoint(bumpIndex);
        bumpSign(nUnit)       = actMat(nUnit, bumpStartPoint(nUnit))/zScores(nUnit, bumpStartPoint(nUnit));
    end
end

[~, similaritySort]       = sortrows([bumpStartPoint, bumpSize, bumpSign], [-3 1 -2]);
% similaritySort            = similaritySort(end:-1:1);
imagesc(zScores(similaritySort,:))
colormap(french(128,2))
caxis([-1 1])

% First separate the dataset into postive/negative z-score




% for nData           = 2:length(fileList)
%     disp(DataSetList(nData).name)
%     load([TempDatDir DataSetList(nData).name '.mat']);
%     ephysCellIndex  = [DataSetList(nData).cellinfo(:).AP_axis]  > 2100 &...
%                       [DataSetList(nData).cellinfo(:).AP_axis]  < 2900 &...
%                       [DataSetList(nData).cellinfo(:).ML_axis]  > 1100 &...
%                       [DataSetList(nData).cellinfo(:).ML_axis]  < 1900;  
%     otherCellIndex  = ~ephysCellIndex; 
% end
