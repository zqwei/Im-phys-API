%
% plotHistActivityExampleUnits.m
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

function barHistWithDist (barData, Dist, xlabels, barSeries, barColor, barSign)

    % Dist  = 'Normal';
    % Dist  = 'Poisson'
    
%     if nargin  < 4
%         Color   = 'k';
%     end
%     
%     if nargin  < 3
%         barSeries = 1000; % numBar in this case
%     end
    
    if nargin<6
        barSign = 1;
    end

    [histFreq, histXout] = hist(barData, barSeries);
    bar(histXout, histFreq*barSign, 'faceColor', barColor, 'edgeColor', barColor);  
    if strcmp(Dist, 'Poisson')
        PD                   = fitdist(floor(barData(barData<10)), Dist); 
    else
        PD                   = fitdist(barData, Dist); 
    end
    if strcmp(Dist, 'Poisson')
        %binWidth             = round(histXout(2)-histXout(1));
        binWidth             = histXout(2)-histXout(1);
        % plot(histXout, pdf(PD, round(histXout))*binWidth*sum(histFreq)*barSign, 'Color', barColor, 'linewid',1.0);
        plot(histXout, pdf(PD, histXout)*binWidth*sum(histFreq)*barSign, 'Color', barColor, 'linewid',1.0);
    else
        binWidth             = histXout(2)-histXout(1);
        plot(histXout, pdf(PD, histXout)*binWidth*sum(histFreq)*barSign, 'Color', barColor, 'linewid',1.0);
    end
    
    ylabel('Count (#)')
    xlabel(xlabels)