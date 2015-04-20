%
%
%
% Ziqiang Wei
% 

function plot3DGPFA(seq, params)

    dimsToPlot = 1:3;
    nPlotMax   = 20;
    randTrace  = randperm(length(seq));
    seq        = seq(randTrace);
    
    yes_trial  = [seq(:).trialType];
    
    [xDim, nT] = size(seq(1).xorth);
    
    xorths     = zeros(length(seq), xDim, nT);

    for  n     = 1:length(seq)
      xorths(n, :, :)  = seq(n).xorth;
    end
    
    timePoints = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    
    figure;
    hold on;
    for n      = 1:nPlotMax
        dat    = squeeze(xorths(n, :, :));
        if yes_trial(n) == 1
            plot3(dat(1,:), dat(2,:), dat(3,:), '-', 'linewidth', 0.5, 'color', [1 0.5 0.5]);
        else
            plot3(dat(1,:), dat(2,:), dat(3,:), '-', 'linewidth', 0.5, 'color', [0.5 0.5 1]);
        end
        
        for nT     = 2:4
            plot3(dat(1,timePoints(nT)), dat(2,timePoints(nT)), dat(3,timePoints(nT)), '.', 'markersize', 10, 'color', [1-0.25*nT 1-0.25*nT 1-0.25*nT]);
        end        
    end
    
    dat    = squeeze(mean(xorths(yes_trial==1, :, :)));
    plot3(dat(1,:), dat(2,:), dat(3,:), '-', 'linewidth', 2, 'color', [1 0 0]);
    for nT     = 2:4
        plot3(dat(1,timePoints(nT)), dat(2,timePoints(nT)), dat(3,timePoints(nT)), '.', 'markersize', 10, 'color', [1-0.25*nT 1-0.25*nT 1-0.25*nT]);
    end
    
    dat    = squeeze(mean(xorths(yes_trial==0, :, :)));
    plot3(dat(1,:), dat(2,:), dat(3,:), '-', 'linewidth', 2, 'color', [0 0 1]);
    for nT     = 2:4
        plot3(dat(1,timePoints(nT)), dat(2,timePoints(nT)), dat(3,timePoints(nT)), '.', 'markersize', 10, 'color', [1-0.25*nT 1-0.25*nT 1-0.25*nT]);
    end
    
    str1 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(1));
    str2 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(2));
    str3 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(3));

    xlabel(str1, 'interpreter', 'latex');
    ylabel(str2, 'interpreter', 'latex');
    zlabel(str3, 'interpreter', 'latex');
    
    view(3)