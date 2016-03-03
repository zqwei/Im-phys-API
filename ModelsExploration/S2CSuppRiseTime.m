%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% low firing activity to calcium analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDirV1Cells

if ~exist([PlotDir 'ModelCellFits'],'dir')
    mkdir([PlotDir 'ModelCellFits'])
end

Fm              = 21.3240;
K               = 13.9248;
n               = 1.7531;
tau_rise        = 0.0728;
tau_decay       = 1.4551;

timeSeriesData  = 0:0.01:4;
numTrial        = 200;
numPlotTrial    = 20;
intNoise        = 1.5;
extNoise        = 1.5;

colorSet        = cbrewer('div', 'Spectral', 20, 'cubic');

frThres         = 2;


for rMean       = 1:20
    Ca              = zeros(numTrial, length(timeSeriesData));
    preSamplePoints = exprnd(1/rMean, numTrial, 100);
    preSamplePoints = cumsum(preSamplePoints, 2);
    
    for nTrial  = 1:numTrial
        for nSpike = 1:100
            Delta_t = max(timeSeriesData - preSamplePoints(nTrial, nSpike),0);
            Ca(nTrial,:) = Ca(nTrial,:) + exp(-Delta_t/tau_decay).*(1-exp(-Delta_t/tau_rise));
        end
    end
    
    intNoiseCa = Ca + randn(size(Ca))*intNoise; 
    intNoiseCa(intNoiseCa<0) = 0;    
    CaImaging  = Fm* Ca.^n ./ (K^n + Ca.^n);
    intNoiseCaImaging = Fm* intNoiseCa.^n ./ (K^n + intNoiseCa.^n);
    extNoiseCaImaging = CaImaging + randn(size(CaImaging))*extNoise; 
    

    subplot(1, 4, 1)
    hold on;
    for nTrial = 1:numPlotTrial
        plot(preSamplePoints(nTrial, :), ones(1, 100)*nTrial + rMean*(numPlotTrial+5), ...
            '.', 'color', colorSet(rMean, :))
    end

    text(0, rMean*(numPlotTrial+5)+13, [num2str(rMean) ' Hz'])

    xlim([timeSeriesData(1) timeSeriesData(end)])
    ylim([numPlotTrial+5 rMean*(numPlotTrial+5)+numPlotTrial+5+0.5])
    axis off
    
    if mod(rMean, 3) == 2
        subplot(1, 4, 2)
        hold on
        plot(timeSeriesData, z1(mean(CaImaging, 1)), ...
            '-', 'color', colorSet(rMean, :), 'linewid', 1);
        box off
        ylabel('normalized F')
        xlabel('Time (sec)')
        ylim([0 1])
        xlim([0 4])
        set(gca, 'TickDir', 'out')
        title ('no noise')

        subplot(1, 4, 3)
        hold on
        plot(timeSeriesData, z1(mean(intNoiseCaImaging, 1)), ...
            '-', 'color', colorSet(rMean, :), 'linewid', 1);
        box off
        ylabel('normalized F')
        xlabel('Time (sec)')
        ylim([0 1])  
        xlim([0 4])
        set(gca, 'TickDir', 'out')
        title ('internal noise')

        subplot(1, 4, 4)
        hold on
        plot(timeSeriesData, z1(mean(extNoiseCaImaging, 1)), ...
            '-', 'color', colorSet(rMean, :), 'linewid', 1);
        box off
        ylabel('normalized F')
        xlabel('Time (sec)')
        ylim([0 1])
        xlim([0 4])
        set(gca, 'TickDir', 'out')
        title ('external noise')
    end
    
end

setPrint(8*4, 8, [PlotDir 'ModelCellFits/S2CSuppRiseTime'])
close all