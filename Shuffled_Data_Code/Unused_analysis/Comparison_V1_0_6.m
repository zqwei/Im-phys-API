% 
% Comparison based on single unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0.5
% 
% Summary for all dataset

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
% Gaussian filter for spiking data
sigma                         = 0.1 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 30;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);

if ~exist([PlotDir '/All_Data_Summary'],'dir')
    mkdir([PlotDir '/All_Data_Summary'])
end
figure;

for nData                     = 1:length(fileList)
    load([TempDatDir DataSetList(nData).name '.mat']);
    numUnits                  = length(nDataSet);
    pValue                    = applyFuncToCompareTrialType(nDataSet, @pValueTTest2);
    meanDiffValue             = applyFuncToCompareTrialType(nDataSet, @meanDiff);
    logPValue                 = -log(pValue);
    % yes  -- blue trial
    % no   -- red trial
    zScores                   = -sign(meanDiffValue).*logPValue;     
    actMat                    = logPValue;
    numT                      = size(zScores,2);
    bumpActThres              = 3; % > bumpActThres considering as a bump % 3 = -log(0.05)
    bumpMat                   = actMat > bumpActThres;
    bumpSize                  = ones(size(actMat,1),1);
    bumpStartPoint            = nan(size(actMat,1),1);
    bumpSign                  = ones(size(actMat,1),1);

    for nUnit = 1: size(actMat,1)
        diffActMat                = diff([bumpMat(nUnit,:),0]);
        beginPoint                = find(diffActMat==1);
        endPoint                  = find(diffActMat==-1);
        if length(endPoint)>length(beginPoint); endPoint = endPoint(2:end); end
        [bumplength, bumpIndex]   = max(endPoint - beginPoint);
        if isempty(bumplength); bumplength = 0; end
        bumpSize(nUnit)           = bumplength;
        if bumplength==0
            bumpStartPoint(nUnit) = numT;
            bumpSign(nUnit)       = -1;
        else 
            bumpStartPoint(nUnit) = beginPoint(bumpIndex);
            bumpSign(nUnit)       = actMat(nUnit, bumpStartPoint(nUnit))/zScores(nUnit, bumpStartPoint(nUnit));
        end
    end
    [~, similaritySort]           = sortrows([bumpStartPoint, bumpSize, bumpSign], [-3 -1 -2]);    

    h = subplot(length(fileList), 4, 4*(nData -1) + 1);
    hold on;
    imagesc(DataSetList(nData).params.timeSeries, 1:numUnits, zScores(similaritySort,:));
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)])
    ylim([1 numUnits])
    colormap(h, french(128,2))
    colorbar;
    caxis([-5 5])
    axis xy
    xlabel('Time (s)')
    ylabel('Neuron Index')
    box off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dependence analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pValue                    = applyFuncToCompareTrialType(nDataSet, @pValueTTest2, filterInUse);
    meanDiffValue             = applyFuncToCompareTrialType(nDataSet, @meanDiff, filterInUse);
    logPValue                 = -log(pValue);
    % yes  -- blue trial
    % no   -- red trial
    zScores                   = -sign(meanDiffValue).*logPValue;     
    actMat                    = logPValue;
    numT                      = size(zScores,2);
    bumpActThres              = 3; % > bumpActThres considering as a bump % 3 = -log(0.05)
    bumpMat                   = actMat > bumpActThres;
    bumpSize                  = ones(size(actMat,1),1);
    bumpStartPoint            = nan(size(actMat,1),1);
    bumpSign                  = ones(size(actMat,1),1);

    for nUnit = 1: size(actMat,1)
        diffActMat                = diff([bumpMat(nUnit,:),0]);
        beginPoint                = find(diffActMat==1);
        endPoint                  = find(diffActMat==-1);
        if length(endPoint)>length(beginPoint); endPoint = endPoint(2:end); end
        [bumplength, bumpIndex]   = max(endPoint - beginPoint);
        if isempty(bumplength); bumplength = 0; end
        bumpSize(nUnit)           = bumplength;
        if bumplength==0
            bumpStartPoint(nUnit) = numT;
            bumpSign(nUnit)       = -1;
        else 
            bumpStartPoint(nUnit) = beginPoint(bumpIndex);
            bumpSign(nUnit)       = actMat(nUnit, bumpStartPoint(nUnit))/zScores(nUnit, bumpStartPoint(nUnit));
        end
    end

%     bumpStartPoint(bumplength<15) = numT;
    ML_axis                      = [DataSetList(nData).cellinfo(:).ML_axis];
    AP_axis                      = [DataSetList(nData).cellinfo(:).AP_axis];
    depth                        = [DataSetList(nData).cellinfo(:).depth];    
    if nData    == 1; depth      = floor((depth-110)/40)*40+110; end
    uniqueDepth                  = unique(depth);
    depthStrings                 = cell(length(uniqueDepth),1);  
    depthStrings(1:end)          = {''};
    if length(uniqueDepth)       <=3
        depthStrings             = cellstr(num2str(uniqueDepth'));
    else
        stepLength               = floor(length(uniqueDepth)/3);
        depthStrings(1:stepLength:end) = cellstr(num2str(uniqueDepth(1:stepLength:end)'));
    end
    
    subplot(length(fileList), 4, 4*(nData -1) + 2)
    hold on;
    iGroup                       = floor((ML_axis-1100)/200);
    numCellInGroup               = histcounts(iGroup, 0:4);
    startTimePoint               = DataSetList(nData).params.timeSeries(bumpStartPoint);
    for nGroup                   = 0:3
        nCell                    = histcounts(startTimePoint(iGroup == nGroup), DataSetList(nData).params.timeSeries(1:3:end-3));
        plot(DataSetList(nData).params.timeSeries(4:3:end-3), cumsum(nCell)/sum(nCell));
        legendString(nGroup + 1) = {['ML ' num2str(1.2+nGroup*0.2) 'mm, ' num2str(numCellInGroup(nGroup+1)) ' cells']}; %#ok<SAGROW>
    end
    legend(legendString, 'Location', 'northwest');
    legend('boxoff')
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)])
    ylim([0 1])
    ylabel('Fraction of Cell')
    xlabel('First Sig. p Val. Time (s)')
    title(strrep(DataSetList(nData).name(9:end),'_',' '))
    
    h = subplot(length(fileList), 4, 4*(nData -1) + 3);
    if nData == 1
        ML_axis                   = ML_axis + (rand(size(ML_axis))-0.5) * 200;
        AP_axis                   = AP_axis + (rand(size(AP_axis))-0.5) * 50;
    end
    
    scatter(ML_axis(bumpStartPoint<numT), AP_axis(bumpStartPoint<numT), [],round(startTimePoint(bumpStartPoint<numT)*3)/3, '.');
    colormap(h, parula)
    colorbar;
    xlim([1000 2100])
    ylim([2100 3000])
    ylabel('AP location (um)')
    xlabel('ML location (um)')
    box off
    
    subplot(length(fileList), 4, 4*(nData -1) + 4)
    if nData == 1
        ML_axis                   = ML_axis + (rand(size(ML_axis))-0.5) * 200;
        AP_axis                   = AP_axis + (rand(size(AP_axis))-0.5) * 50;
    end
    
    boxplot(startTimePoint(bumpStartPoint<numT), depth(bumpStartPoint<numT),'labels',depthStrings);
%     xlim([110 750])
%     ylim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)])
    ylabel('First Sig. p Val. Time (s)')
    xlabel('Depth (um)')
    box off
end

setPrint(4*16, length(DataSetList)*12, [PlotDir 'All_Data_Summary/All_Data_Summary'], 'pdf')