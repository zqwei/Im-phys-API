%
% plotAnovaLogP.m
%
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

function plotAnovaLogP (logPValue, params)

    figure;
    thresLogP     = -log(0.05);
    numTime       = 10;
    numUnit       = size(logPValue, 1);
    
    selectivityMat = sum(logPValue > thresLogP, 3) > numTime;
    unitGroup     = zeros(numUnit, 1);
    
    gAP           = selectivityMat(:, 1);
    gLR           = selectivityMat(:, 2);
    gCE           = selectivityMat(:, 3);
    
    unitGroup( gAP & ~gLR & ~gCE) = 1; % group for AP
    unitGroup(~gAP &  gLR & ~gCE) = 2; % group for LR
    unitGroup(~gAP & ~gLR &  gCE) = 3; % group for CE
    unitGroup( gAP &  gLR & ~gCE) = 4; % group for AP + LR
    unitGroup( gAP & ~gLR &  gCE) = 5; % group for AP + CE
    unitGroup(~gAP &  gLR &  gCE) = 6; % group for LR + CE
    unitGroup( gAP &  gLR &  gCE) = 7; % group for AP + LR + CE 
    unitGroup(~gAP & ~gLR & ~gCE) = 8; % group for none
    
    [~, unitIndex] = sort(unitGroup, 'ascend');
    sizeGroup      = histcounts(unitGroup, 1:9);
    groupTitle     = {'-log(P) for pole location', '-log(P) for lick direction', '-log(P) for reward'};
    
    for nGroup     = 1:3
        subplot(1, 4, nGroup)
        hold on;
        imagesc(params.timeSeries, 1:numUnit, squeeze(logPValue(unitIndex, nGroup, :)));
        caxis([0 5])
        gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
        hold off;
        xlim([params.timeSeries(1) params.timeSeries(end)]);
        ylim([1 numUnit])
        xlabel('Time (s)');
        ylabel('Neuronal index')
        title(groupTitle{nGroup});
    end
    
    subplot(1, 4, 4)
    pie(sizeGroup, {'pole', 'lick', 'reward', 'pole+lick', 'pole+reward', 'lick+reward', 'all', 'none'})
    
end

