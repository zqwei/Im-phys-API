%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

ylabels                 = {'Fring Rate (Hz)', 'dF/F', 'dF/F', 'dF/F', 'dR/R' };
yAxes_set               = [0 60; -0.5 2.0; -0.5 2.0; -0.5 2.0; -0.5 2.0 ; -0.5 2.0];
lowFiringThres          = [15, 0.3, 0.3, 0.3, 0.3, 0.3];

if ~exist([PlotDir 'ModeledComparingSingleUnitsImagescWithSort'],'dir')
    mkdir([PlotDir 'ModeledComparingSingleUnitsImagescWithSort'])
end


nData = 1;
load([TempDatDir DataSetList(nData).name '.mat'])

blankSpace = 10;
numT       = size(nDataSet(1).unit_yes_trial,2);
actMat     = nan(length(nDataSet),numT*2+blankSpace);

for nUnit = 1: length(nDataSet)
    actMat(nUnit, 1:numT)     = mean(nDataSet(nUnit).unit_yes_trial);
    actMat(nUnit, numT+blankSpace+1:end) = mean(nDataSet(nUnit).unit_no_trial);
end

maxActMat                     = nanmax(actMat, [], 2);
minActMat                     = nanmin(actMat, [], 2);

actMat                    = bsxfun(@minus, actMat, minActMat);
tMinActMat                = minActMat;
actMat                    = bsxfun(@rdivide, actMat, maxActMat-tMinActMat);


bumpActThres              = 0.6;
bumpMat                   = actMat > bumpActThres;
bumpSize                  = ones(size(actMat,1),1);
bumpStartPoint            = nan(size(actMat,1),1);

for nUnit = 1: size(actMat,1)
    diffActMat            = diff([bumpMat(nUnit,:),0]);
    beginPoint            = find(diffActMat==1);
    endPoint              = find(diffActMat==-1);
    if length(endPoint)>length(beginPoint); endPoint = endPoint(2:end); end
    [bumplength, bumpIndex] = max(endPoint - beginPoint);
    if isempty(bumplength); bumplength = 0; end
    bumpSize(nUnit)       = bumplength;
    if bumplength==0
        bumpStartPoint(nUnit) = numT*2+blankSpace;
    else
        bumpStartPoint(nUnit) = beginPoint(bumpIndex);
    end
end

[~, similaritySort]       = sortrows([bumpStartPoint, bumpSize], [1 -2]);
similaritySort            = similaritySort(end:-1:1);
% 
% figure;
% hold on;
% timeSeries                    =  DataSetList(nData).params.timeSeries;
% minTime                       =  DataSetList(nData).params.timeSeries(1);
% maxTime                       =  DataSetList(nData).params.timeSeries(end);
% betweenSpace                  = blankSpace * DataSetList(nData).params.binsize;
% constShift                    =  - minTime + maxTime + betweenSpace;
% timeSeries                    = [timeSeries, maxTime + (1:blankSpace)*DataSetList(nData).params.binsize, timeSeries+constShift];
% tMaxTime                      = timeSeries(end);
% b = imagesc(timeSeries, 1:length(similaritySort), actMat(similaritySort, :), [0 1]);
% set(b,'AlphaData',~isnan(actMat));
% axis xy;
% %gridxy ([maxTime + betweenSpace/2],[], 'Color','k','Linestyle','-','linewid', 2.0)
% gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0, DataSetList(nData).params.polein+constShift, DataSetList(nData).params.poleout+constShift, constShift],[], 'Color','k','Linestyle','--','linewid', 1.0)
% ax              = gca;
% hColor          = colorbar;    
% ylabel(hColor, 'Normalized activity');
% cPos            = hColor.Position;
% axpos           = ax.Position;
% hColor.Position = [axpos(1)+axpos(3)+0.11 cPos(2)+cPos(4)*0.25 cPos(3)*0.5 cPos(4)*0.5];    
% box off
% xTicks                        = round(minTime):floor(maxTime);
% xTickLabel                    = arrayfun(@(x) num2str(x), xTicks,'Uniform',false);
% xTickLabel(mod(xTicks,2)==1)  = {''};
% set(ax, 'XTick', [xTicks xTicks+constShift], 'XTickLabel', [xTickLabel, xTickLabel]);
% set(ax, 'YTick', [1 length(nDataSet)])
% axis([minTime, tMaxTime, 1, length(similaritySort)]);
% 
% xlabel('Time (s)')
% ylabel('Neuronal index');
% 
% hold off
% 
% setPrint(6*2, 3*3, [PlotDir 'ModeledComparingSingleUnitsImagescWithSort/SingleUnitsImagescWithSort_' DataSetList(nData).name])





load ([TempDatDir 'DataListS2CModel.mat']);
for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    
    for nUnit = 1: length(nDataSet)
        actMat(nUnit, 1:numT)     = mean(nDataSet(nUnit).unit_yes_trial);
        actMat(nUnit, numT+blankSpace+1:end) = mean(nDataSet(nUnit).unit_no_trial);
    end

    maxActMat                     = nanmax(actMat, [], 2);
    minActMat                     = nanmin(actMat, [], 2);

    actMat                    = bsxfun(@minus, actMat, minActMat);
    tMinActMat                = minActMat;
    actMat                    = bsxfun(@rdivide, actMat, maxActMat-tMinActMat);
    
    
    figure;
    hold on;
    b = imagesc(timeSeries, 1:length(similaritySort), actMat(similaritySort, :), [0 1]);
    set(b,'AlphaData',~isnan(actMat));
    axis xy;
    %gridxy ([maxTime + betweenSpace/2],[], 'Color','k','Linestyle','-','linewid', 2.0)
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0, DataSetList(nData).params.polein+constShift, DataSetList(nData).params.poleout+constShift, constShift],[], 'Color','k','Linestyle','--','linewid', 1.0)
    ax              = gca;
    hColor          = colorbar;    
    ylabel(hColor, 'Normalized activity');
    cPos            = hColor.Position;
    axpos           = ax.Position;
    hColor.Position = [axpos(1)+axpos(3)+0.11 cPos(2)+cPos(4)*0.25 cPos(3)*0.5 cPos(4)*0.5];    
    box off
    xTicks                        = round(minTime):floor(maxTime);
    xTickLabel                    = arrayfun(@(x) num2str(x), xTicks,'Uniform',false);
    xTickLabel(mod(xTicks,2)==1)  = {''};
    set(ax, 'XTick', [xTicks xTicks+constShift], 'XTickLabel', [xTickLabel, xTickLabel]);
    set(ax, 'YTick', [1 length(nDataSet)])
    axis([minTime, tMaxTime, 1, length(similaritySort)]);

    xlabel('Time (s)')
    ylabel('Neuronal index');

    hold off

    setPrint(6*2, 3*3, [PlotDir 'ModeledComparingSingleUnitsImagescWithSort/SingleUnitsImagescWithSort_' DataSetList(nData).name])
        
end

close all;