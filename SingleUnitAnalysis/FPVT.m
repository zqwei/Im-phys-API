%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first p-value ramping time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
% Gaussian filter for spiking data
sigma                         = 0.1 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 10;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);

if ~exist([PlotDir 'SingleUnitsFPVT'],'dir')
    mkdir([PlotDir 'SingleUnitsFPVT'])
end


for nData                     = 1:length(DataSetList)-1
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

    h = figure;
    hold on
    imagesc(DataSetList(nData).params.timeSeries, 1:numUnits, zScores(similaritySort,:));
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)])
    ylim([1 numUnits])
    colormap(h, french(128,2))
%     colorbar;
    caxis([-5 5])
    axis xy
    xlabel('Time (s)')
    ylabel('Neuron Index')
    box off;
    setPrint(8, 6, [PlotDir 'SingleUnitsFPVT/SingleUnitsZScore_' DataSetList(nData).name], 'pdf')

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

    % bumpStartPoint(bumplength<15) = numT;
    ML_axis                      = [DataSetList(nData).cellinfo(:).ML_axis];
    AP_axis                      = [DataSetList(nData).cellinfo(:).AP_axis];
    depth                        = [DataSetList(nData).cellinfo(:).depth];
    
    depthStart                   = 100;
    depthBin                     = 50;
    depthEnd                     = 900;    
    depth                        = floor((depth-depthStart)/depthBin)*depthBin+depthStart;
    depth(depth>depthEnd)        = depthEnd;
    depth(depth<depthStart)      = depthStart;
    
%     uniqueDepth                  = unique(depth);
    uniqueDepth                  = depthStart:depthBin:depthEnd;
    depthStrings                 = cell(length(uniqueDepth),1);  
    depthStrings(1:end)          = {''};
    if length(uniqueDepth)       <=3
        depthStrings             = cellstr(num2str(uniqueDepth'));
    else
        stepLength               = floor(length(uniqueDepth)/3);
        depthStrings(1:stepLength:end) = cellstr(num2str(uniqueDepth(1:stepLength:end)'));
    end

    iGroup                       = floor((ML_axis-1100)/200);
    numCellInGroup               = histcounts(iGroup, 0:4);
    startTimePoint               = DataSetList(nData).params.timeSeries(bumpStartPoint);

    
    
    
    timeBins                     = DataSetList(nData).params.timeSeries;
    pdfMat                       = nan(length(timeBins), length(uniqueDepth));
    numCells                     = zeros(length(uniqueDepth),1);
    
    
    for nDepth = 1:length(uniqueDepth)
        numCells(nDepth)         = sum(bumpStartPoint<numT & depth' == uniqueDepth(nDepth));
        if sum(bumpStartPoint<numT & depth' == uniqueDepth(nDepth))>0
            pdfMat(:, nDepth)        = ksdensity(startTimePoint(bumpStartPoint<numT & depth' == uniqueDepth(nDepth)), timeBins);
        end
    end
    
    figure;
    subplot(1, 2, 1)
    hold on
    h = imagesc(timeBins, uniqueDepth, pdfMat');
    set(h,'alphadata',~isnan(pdfMat'))
    axis xy
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)])
    ylim([0 950])
    set(gca, 'yTick', 0:300:900)
    colorbar
    xlabel('First Sig. p Val. Time (s)')
    ylabel('Depth (um)')
    title('Empricial prob. density')
    box off
    
    subplot(1, 2, 2)
    bar(uniqueDepth, numCells)
    ylabel('# cells')
    xlabel('Depth (um)')
    xlim([0 950])
    set(gca, 'xTick', 0:300:900)
    
    
    setPrint(8*2, 6, [PlotDir 'SingleUnitsFPVT/SingleUnitsFPVTDepth_' DataSetList(nData).name], 'pdf')
end

close all