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

function plotMeanActivityImagescWithSort (nDataSet, params, maxValue, minValue, ylabels, lowFiringThres, yAxes_set)    
    hDendrogram = figure;
    blankSpace = 10;
    numT       = size(nDataSet(1).unit_yes_trial,2);
    actMat     = nan(length(nDataSet),numT*2+blankSpace);
    
    for nUnit = 1: length(nDataSet)
        actMat(nUnit, 1:numT)     = mean(nDataSet(nUnit).unit_yes_trial);
        actMat(nUnit, numT+blankSpace+1:end) = mean(nDataSet(nUnit).unit_no_trial);
    end
    
    % 1. normalization of the neuronal activity to range [0, 1]
    
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
    
    % 2. sort the normalized activity using linkage; parameters 'single', 'corrleation'
    similarityValue               = linkage(actMat(:, [1:numT, numT+blankSpace+1:end]), 'complete','correlation');
    [~, ~, similaritySort]        = dendrogram(similarityValue,0);
%     D                             = pdist(actMat(:, [1:numT, numT+blankSpace+1:end]),'correlation');
%     similaritySort                = optimalleaforder(similarityValue,D, 'Criteria', 'group');

    similaritySort                = similaritySort(end:-1:1);
    

    close(hDendrogram);
    figure;
    % 3. plot of imagesc
    subplot(3, 6, [1 2 7 8 13 14])
    hold on;
    timeSeries                    = params.timeSeries;
    minTime                       = params.timeSeries(1);
    maxTime                       = params.timeSeries(end);
    betweenSpace                  = blankSpace * params.binsize;
    constShift                    =  - minTime + maxTime + betweenSpace;
    timeSeries                    = [timeSeries, maxTime + (1:blankSpace)*params.binsize, timeSeries+constShift];
    tMaxTime                      = timeSeries(end);
    b = imagesc(timeSeries, 1:length(nDataSet), actMat(similaritySort, :), [0 1]);
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
    axis([minTime, tMaxTime, 1, length(nDataSet)]);
    
    xlabel('Time (s)')
    ylabel('Neuronal index');
    
    hold off
    
    % 4. plot of example neurons
    numExampleNeurons            = 9;
    % only using the neurons with a high activity
    %
    % numExampleNeurons is also the number of the clusters
    %
    % 4.1 filtering out the low firing units
    % lowFiringThres                = 0.5;
    highActMatIndex               = find(maxActMat> lowFiringThres);
    highActMat                    = actMat(highActMatIndex, :);
%     clusterIndex                 = clusterWithKickOut(highActMat(:, [1:numT, numT+blankSpace+1:end]), numExampleNeurons, 2);
    clusterIndex                  = cluster(...
                                linkage(highActMat(:, [1:numT, numT+blankSpace+1:end]), 'single','correlation'),...
                                'maxclust',numExampleNeurons);

%     m                             = ceil(numExampleNeurons/3);
%     nCluster                      = 1;
%     for nNeuron                   = 1:numExampleNeurons        
%         if sum(clusterIndex==nNeuron) > 2 && nCluster <= 9
%             subplot(3, 6, (mod(nCluster-1,3)+4) + floor((nCluster-1)/3)*6)
%             nNeuronIndex              = highActMatIndex(find(clusterIndex==nNeuron,1,'first')); 
%             plotMeanActivityTrace (nDataSet, nNeuronIndex, params, ylabels)
%             nCluster                  = nCluster + 1;
%         end
%     end

%     m                             = ceil(numExampleNeurons/3);

%     timePoints                    = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);  
%     negResponseNeuronIndex        = mean(actMat(:,timePoints(4):timePoints(5)),2)<0;
%     negResponseNeuronIndex        = find(negResponseNeuronIndex);
    
%     figure;
%     
%     m = ceil(sqrt(length(negResponseNeuronIndex)));
%     
%     for nNeuron                   = 1: length(negResponseNeuronIndex)% numExampleNeurons 
%         subplot(m, m, nNeuron)
%         nNeuronIndex              = negResponseNeuronIndex(nNeuron);
%         plotMeanActivityTrace (nDataSet, nNeuronIndex, params, ylabels, ' ', yAxes_set)
%     end

%     for nNeuron                   = 1:numExampleNeurons        
%         subplot(3, 6, (mod(nNeuron-1,3)+4) + floor((nNeuron-1)/3)*6)
%         nNeuronIndex              = negResponseNeuronIndex(nNeuron); % highActMatIndex(find(clusterIndex==nNeuron,1,'first')); 
%         if nNeuron == 4
%             plotMeanActivityTrace (nDataSet, nNeuronIndex, params, ylabels, ' ', yAxes_set)
%         elseif nNeuron == 8
%             plotMeanActivityTrace (nDataSet, nNeuronIndex, params, ' ', 'Time (s)', yAxes_set)
%         else
%             plotMeanActivityTrace (nDataSet, nNeuronIndex, params, ' ', ' ', yAxes_set)
%         end
%     end
    
end

function plotMeanActivityTrace (nDataSet, nNeuron, params, ylabels, xlabels, yAxes_set)
    mean_yes      = mean(nDataSet(nNeuron).unit_yes_trial);
    mean_no       = mean(nDataSet(nNeuron).unit_no_trial);

    var_yes       = sem(nDataSet(nNeuron).unit_yes_trial);
    var_no        = sem(nDataSet(nNeuron).unit_no_trial);
    

    hold on;
    shadedErrorBar(params.timeSeries, mean_yes, var_yes, {'-b', 'linewid', 1.0}, 0.5);
    shadedErrorBar(params.timeSeries, mean_no, var_no, {'-r', 'linewid', 1.0}, 0.5);
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    % ylim([yAxes_set(1) yAxes_set(2)]);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off; 
    xlabel(xlabels);
    ylabel(ylabels)
    
end
