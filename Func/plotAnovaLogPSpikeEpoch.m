%
% plotAnovaLogPSpike.m
% 
%
% Spiking dataset
%
% ----------------------------
% Output:
%
% ----------------------------
% version 1.0
%
% ----------------------------
% version 1.1
% 
% - get group infomation from getAnovaLogPSpikeEpoch.m
% - remove pie plot
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function unitGroup = plotAnovaLogPSpikeEpoch (logPValue, unitGroup, params)

    figure;
    cmap = cbrewer('qual', 'Set2', 10, 'cubic'); 
    
    [sortedIndex, unitIndex] = sort(unitGroup, 'ascend');
    groupTitle     = {  '-log(P) for pole location', ...
                        '-log(P) for lick direction', ...
                        '-log(P) for reward'};
    
    sortLogPValue  = logPValue(unitIndex, :, :);
    sortLogPValue  = sortLogPValue(sortedIndex < 8, :, :);
    numUnit        = size(sortLogPValue, 1);
    
    for nGroup     = 1:3
        subplot(3, 10, (nGroup-1)*10 + 1)
        imagesc(sortedIndex(sortedIndex<8));
        caxis([1 8])
        axis xy
        axis off
        colormap(cmap)
        freezeColors
        
        subplot(3, 10, (nGroup-1)*10 + (4:10))
        hold on;
        imagesc(params.timeSeries, 1:numUnit, squeeze(sortLogPValue(:, nGroup, :)));
        caxis([0 4])
        gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
        hold off;
        colorbar
        xlim([params.timeSeries(1) params.timeSeries(end)]);
        ylim([1 numUnit])
        xlabel('Time (s)');
        ylabel('Neuronal index')
        set(gca, 'YTick', [1 numUnit])
        set(gca, 'TickDir', 'out')
        title(groupTitle{nGroup});
        colormap gray
        freezeColors

        
    end
    
    
end


% version 1.0 code
%
% function unitGroup = plotAnovaLogPSpikeEpoch (logPValue, logPValueEpoch, params)
% 
%     figure;
%     thresLogP     = -log(0.05);
%     numUnit       = size(logPValue, 1);
%     unitGroup     = nan(numUnit, 1);
% 
%     gAP           = sum(logPValueEpoch(:, 1, :) > thresLogP, 3) > 0;
%     gLR           = sum(logPValueEpoch(:, 2, :) > thresLogP, 3) > 0;
%     gCE           = sum(logPValueEpoch(:, 3, 3) > thresLogP, 3) > 0;
% 
% 
%     unitGroup( gAP & ~gLR & ~gCE) = 1; % group for AP
%     unitGroup(~gAP &  gLR & ~gCE) = 2; % group for LR
%     unitGroup(~gAP & ~gLR &  gCE) = 3; % group for CE
%     unitGroup( gAP &  gLR & ~gCE) = 4; % group for AP + LR
%     unitGroup( gAP & ~gLR &  gCE) = 5; % group for AP + CE
%     unitGroup(~gAP &  gLR &  gCE) = 6; % group for LR + CE
%     unitGroup( gAP &  gLR &  gCE) = 7; % group for AP + LR + CE 
%     unitGroup(~gAP & ~gLR & ~gCE & ~isnan(sum(sum(logPValueEpoch, 3), 2))) = 8; % group for none
% %     unitGroup(sum(logPValueEpoch(:, 1, :) <= thresLogP, 3)==3 & ...
% %               sum(logPValueEpoch(:, 2, :) <= thresLogP, 3)==3 & ...
% %               sum(logPValueEpoch(:, 3, 3) <= thresLogP, 3)==1) = 8;
%     
%     [sortedIndex, unitIndex] = sort(unitGroup, 'ascend');
%     sizeGroup      = histcounts(unitGroup, 1:9);
%     groupTitle     = {  '-log(P) for pole location', ...
%                         '-log(P) for lick direction', ...
%                         '-log(P) for reward'};
%     
%     sortLogPValue  = logPValue(unitIndex, :, :);
%     sortLogPValue  = sortLogPValue(sortedIndex < 8, :, :);
%     numUnit        = size(sortLogPValue, 1);
%     
%     for nGroup     = 1:3
%         subplot(4, 10, (nGroup-1)*10 + 1)
%         imagesc(sortedIndex(sortedIndex<8));
%         caxis([1 8])
%         axis xy
%         axis off
%         subplot(4, 10, (nGroup-1)*10 + (4:10))
%         hold on;
%         imagesc(params.timeSeries, 1:numUnit, squeeze(sortLogPValue(:, nGroup, :)));
%         caxis([0 4])
%         gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
%         hold off;
%         colorbar
%         xlim([params.timeSeries(1) params.timeSeries(end)]);
%         ylim([1 numUnit])
%         xlabel('Time (s)');
%         ylabel('Neuronal index')
%         set(gca, 'YTick', [1 numUnit])
%         set(gca, 'TickDir', 'out')
%         title(groupTitle{nGroup});
%     end
%     
%     subplot(4, 10, 31:40)
%     groupNames      = {'Pole', 'Lick', 'Reward', 'PL', 'PR', 'LR', 'PLR', 'Non-Selective'};
%     pie(sizeGroup, groupNames)
%     title('Distribution of cell types')
%     
% end

