%
% plotAnovaLogPSpike.m
% 
%
% Spiking dataset
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

function unitGroup = plotAnovaLogPSpike (logPValue, params)

    figure;
    thresLogP     = -log(0.05);
    numTime       = 15;
    numUnit       = size(logPValue, 1);
    timePoints    = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
%     selectivityMat = sum(logPValue(:, :, timePoints(2):end) > thresLogP, 3) > numTime;
    unitGroup     = zeros(numUnit, 1);
    
%     gAP           = selectivityMat(:, 1);
%     gLR           = selectivityMat(:, 2);
%     gCE           = selectivityMat(:, 3);

    gAP           = sum(logPValue(:, 1, timePoints(2):timePoints(4)) > thresLogP, 3) > numTime;
    gLR           = sum(logPValue(:, 2, timePoints(2):timePoints(4)) > thresLogP, 3) > numTime;
    gCE           = sum(logPValue(:, 3, timePoints(4):timePoints(5)) > thresLogP, 3) > numTime;


    unitGroup( gAP & ~gLR & ~gCE) = 1; % group for AP
    unitGroup(~gAP &  gLR & ~gCE) = 2; % group for LR
    unitGroup(~gAP & ~gLR &  gCE) = 3; % group for CE
    unitGroup( gAP &  gLR & ~gCE) = 4; % group for AP + LR
    unitGroup( gAP & ~gLR &  gCE) = 5; % group for AP + CE
    unitGroup(~gAP &  gLR &  gCE) = 6; % group for LR + CE
    unitGroup( gAP &  gLR &  gCE) = 7; % group for AP + LR + CE 
    unitGroup(~gAP & ~gLR & ~gCE) = 8; % group for none
    
    [sortedIndex, unitIndex] = sort(unitGroup, 'ascend');
    sizeGroup      = histcounts(unitGroup, 1:9);
    groupTitle     = {  '-log(P) for pole location', ...
                        '-log(P) for lick direction', ...
                        '-log(P) for reward'};
    
    sortLogPValue  = logPValue(unitIndex, :, :);
    sortLogPValue  = sortLogPValue(sortedIndex < 8, :, :);
    numUnit        = size(sortLogPValue, 1);
    
    for nGroup     = 1:3
        subplot(4, 10, (nGroup-1)*10 + 1)
        imagesc(sortedIndex(sortedIndex<8));
        caxis([1 8])
        axis xy
        axis off
        subplot(4, 10, (nGroup-1)*10 + (3:10))
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
        title(groupTitle{nGroup});
    end
    
    subplot(4, 10, 31:40)
    groupNames      = {'Pole', 'Lick', 'Reward', 'PL', 'PR', 'LR', 'PLR', 'Other'};
    pie(sizeGroup, groupNames)
    title('Distribution of cell types')
    
end

