%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numComps       = 10;

if ~exist([PlotDir 'CollectedUnitsPCA'],'dir')
    mkdir([PlotDir 'CollectedUnitsPCA'])
end

for nData              = [1 3 4]
    if nData   == 1
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end
    firingRates        = generateDPCAData(nDataSet, numTrials);
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));
    pcaFiringRatesAverage = zeros(numComps, 2, 77);
    firingRatesAverage = [squeeze(firingRatesAverage(:, 1, :)), squeeze(firingRatesAverage(:, 2, :));];
    [coeff,score,~]        = pca(firingRatesAverage', 'NumComponents', numComps);
    actMatTemp             = coeff * score';

    numT       = size(nDataSet(1).unit_yes_trial,2);
    if numT <= 77
        blankSpace = 10;
    else
        blankSpace = round(10/77*numT);
    end
    actMat         = nan(length(nDataSet),numT*2+blankSpace);
    positivePeak   = true(length(nDataSet),1);
    actMat(:, 1:numT)     = actMatTemp(:, 1:numT);
    actMat(:, numT+blankSpace+1:end) = actMatTemp(:, 1+numT:end);
%     for nUnit = 1: length(nDataSet)        
%         positivePeak(nUnit)       = mean(actMat(1:8)) <= mean(actMat(9:47)) ...
%                                    || mean(actMat((1:8)+numT)) <= mean(actMat((9:47)+numT));
%     end
    actMat                        = actMat(positivePeak, :);
    maxActMat                     = nanmax(actMat, [], 2);
    minActMat                     = nanmin(actMat, [], 2);
    actMat                    = bsxfun(@minus, actMat, minActMat);
    tMinActMat                = minActMat;
    actMat                    = bsxfun(@rdivide, actMat, maxActMat-tMinActMat);

    
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
        
    [~, similaritySort]       = sortrows(bumpStartPoint, 1);%sortrows([bumpStartPoint, bumpSize], [1 -2]);
    similaritySort            = similaritySort(end:-1:1);
    params                    = DataSetList(nData).params;
    
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
    set(b,'AlphaData',~isnan(actMat));
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
    set(ax, 'YTick', [1 size(actMat, 1)])
    axis([minTime, tMaxTime, 1, length(similaritySort)]);
    set(ax, 'TickDir', 'Out')
    xlabel('Time (s)')
    ylabel('Neuronal index');
    
    hold off
    
end

% close all