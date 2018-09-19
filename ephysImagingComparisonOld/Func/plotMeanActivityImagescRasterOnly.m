%
% plotMeanActivityImagescWithSort.m
% 
% based on plotMeanActivityImagesc.m
%
% ----------------------------
% Output:
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function plotMeanActivityImagescRasterOnly (nDataSet, params, maxValue, minValue, ylabels)    

    
    
    numT       = size(nDataSet(1).unit_yes_trial,2);
    if numT <= 77
        blankSpace = 10;
    else
        blankSpace = round(10/77*numT);
    end
    
    
    actMat     = nan(length(nDataSet),numT*2+blankSpace);
    
    for nUnit = 1: length(nDataSet)
        actMat(nUnit, 1:numT)     = mean(nDataSet(nUnit).unit_yes_trial, 1);
        actMat(nUnit, numT+blankSpace+1:end) = mean(nDataSet(nUnit).unit_no_trial, 1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. normalization of the neuronal activity to range [0, 1]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    maxActMat                     = prctile(actMat, 98, 2); % nanmax(actMat, [], 2);
    minActMat                     = prctile(actMat, 2, 2); % nanmin(actMat, [], 2);

%     maxActMat                     = nanmax(actMat, [], 2);
%     minActMat                     = nanmin(actMat, [], 2);
    
    if ~isempty(minValue)
        actMat                    = actMat-minValue;
        tMinActMat                = minValue;
    else
        actMat                    = bsxfun(@minus, actMat, minActMat);
        tMinActMat                = minActMat;
    end
    
    if ~isempty(maxValue)
        actMat                    = bsxfun(@rdivide, actMat, maxValue-tMinActMat);
    else
        actMat                    = bsxfun(@rdivide, actMat, maxActMat-tMinActMat);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. sort the neurons
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    bumpActThres              = 0.6; % > bumpActThres considering as a bump
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
    
    figure;
    % 3. plot of imagesc
    hold on;
    timeSeries                    = params.timeSeries;
    minTime                       = params.timeSeries(1);
    maxTime                       = params.timeSeries(end);
    betweenSpace                  = blankSpace / params.frameRate;
    constShift                    =  - minTime + maxTime + betweenSpace;
    timeSeries                    = [timeSeries, maxTime + (1:blankSpace)*params.binsize, timeSeries+constShift];
    tMaxTime                      = timeSeries(end);
    b = imagesc(timeSeries, 1:length(similaritySort), actMat(similaritySort, :), [0 1]);
    set(b,'AlphaData',~isnan(actMat(similaritySort, :)));
    axis xy;
    %gridxy ([maxTime + betweenSpace/2],[], 'Color','k','Linestyle','-','linewid', 2.0)
    gridxy ([params.polein, params.poleout, 0, params.polein+constShift, params.poleout+constShift, constShift],[], 'Color','k','Linestyle','--','linewid', 1.0)
    ax              = gca;
    hColor          = colorbar;    
    ylabel(hColor, 'Normalized activity');
    cPos            = hColor.Position;
    axpos           = ax.Position;
%     ax.Position     = axpos;
    hColor.Position = [axpos(1)+axpos(3)+0.11 cPos(2)+cPos(4)*0.25 cPos(3)*0.5 cPos(4)*0.5];    
    box off
    xTicks                        = round(minTime):floor(maxTime);
    xTickLabel                    = arrayfun(@(x) num2str(x), xTicks,'Uniform',false);
    xTickLabel(mod(xTicks,2)==1)  = {''};
    set(ax, 'XTick', [xTicks xTicks+constShift], 'XTickLabel', [xTickLabel, xTickLabel]);
    set(ax, 'YTick', [1 length(nDataSet)])
    axis([minTime, tMaxTime, 1, length(similaritySort)]);
    set(ax, 'TickDir', 'Out')
    xlabel('Time (s)')
    ylabel('Neuronal index');
    
    hold off
    
end
