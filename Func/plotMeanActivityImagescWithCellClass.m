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

function plotMeanActivityImagescWithCellClass(nDataSet, params, cellinfo, minValue, maxValue)    

    blankSpace = 10;
    numT       = size(nDataSet(1).unit_yes_trial,2);
    numUnit    = length(nDataSet);
    actMat     = nan(numUnit,numT*2+blankSpace);
    
    for nUnit = 1: numUnit
        actMat(nUnit, 1:numT)                = mean(nDataSet(nUnit).unit_yes_trial);
        actMat(nUnit, numT+blankSpace+1:end) = mean(nDataSet(nUnit).unit_no_trial);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. normalization of the neuronal activity to range [0, 1]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    maxActMat                     = nanmax(actMat, [], 2);
    minActMat                     = nanmin(actMat, [], 2);

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
    % Only analyzing the prymidal cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    actMat                        = actMat([cellinfo.cellType]==1,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. classify the neurons
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Separate the tasks into 4 periods
    %
    timePoints                    = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    % 
    %
    % Pole location cells -- classifty pole location regardless decision (
    % correct / error)
    %
    % t-test evaluation (since only two levels: posterior and anterior)
    %
    % 
    locMat                        = nan(numUnit, length(timePoints)-1);
    for nUnit = 1:numUnit
        for nPeriod = 1:length(timePoints)-1
            locMat(nUnit, nPeriod)= ttest2(mean([nDataSet(nUnit).unit_yes_trial(:,timePoints(nPeriod):timePoints(nPeriod+1));...
                                           nDataSet(nUnit).unit_yes_error(:,timePoints(nPeriod):timePoints(nPeriod+1))],2),...
                                           mean([nDataSet(nUnit).unit_no_trial(:,timePoints(nPeriod):timePoints(nPeriod+1));...
                                           nDataSet(nUnit).unit_no_error(:,timePoints(nPeriod):timePoints(nPeriod+1))],2));
        end
    end
    
    locMat(isnan(locMat))         = 0;
    % 
    %
    % Decision cells -- classifty decision regardless pole position (
    % correct / error)
    %
    % t-test evaluation (since only two levels: left and right)
    %
    % 
    lickMat                       = nan(numUnit, length(timePoints)-1);
    for nUnit = 1:numUnit
        for nPeriod = 1:length(timePoints)-1
            lickMat(nUnit, nPeriod)= ttest2(mean([nDataSet(nUnit).unit_yes_trial(:,timePoints(nPeriod):timePoints(nPeriod+1));...
                                           nDataSet(nUnit).unit_no_error(:,timePoints(nPeriod):timePoints(nPeriod+1))],2),...
                                           mean([nDataSet(nUnit).unit_no_trial(:,timePoints(nPeriod):timePoints(nPeriod+1));...
                                           nDataSet(nUnit).unit_yes_error(:,timePoints(nPeriod):timePoints(nPeriod+1))],2));
        end
    end
    
    lickMat(isnan(lickMat))        = 0;
    

    figure;
    % 3. plot of imagesc
    hold on;
    timeSeries                    = params.timeSeries;
    minTime                       = params.timeSeries(1);
    maxTime                       = params.timeSeries(end);
    betweenSpace                  = blankSpace * params.binsize;
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
    axis([minTime, tMaxTime, 1, length(similaritySort)]);
    
    xlabel('Time (s)')
    ylabel('Neuronal index');
    
    hold off
    
    
end
