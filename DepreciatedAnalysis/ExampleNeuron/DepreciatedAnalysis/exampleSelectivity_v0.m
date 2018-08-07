%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% low firing activity to calcium analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir

if ~exist([PlotDir 'ModelCellFits'],'dir')
    mkdir([PlotDir 'ModelCellFits'])
end

Fm              = 1;%21.3240;
Ca0             = 4;
beta            = 1.5;
tau_rise        = 0.0728;
tau_decay       = 3.4551;

timeSeriesData  = 0:0.01:4;
numTrial        = 200;
numPlotTrial    = 20;
intNoise        = 0.15;
extNoise        = 1.5;

colorSet        = cbrewer('div', 'Spectral', 10, 'cubic');

frThres         = 2;

% figure;
for rMean       = 1:10
    Ca              = zeros(numTrial, length(timeSeriesData));
    preSamplePoints = exprnd(1/rMean, numTrial, 100);
    preSamplePoints = cumsum(preSamplePoints, 2);
    preSamplePoints = preSamplePoints(preSamplePoints(nTrial, nSpike) < 0.5);
    maxPre          = max(preSamplePoints);
    preSamplePoints2 = exprnd(1, numTrial, 100);
    
    
    for nTrial  = 1:numTrial
        for nSpike = 1:100
            if preSamplePoints(nTrial, nSpike) < 0.5
                Delta_t = max(timeSeriesData - preSamplePoints(nTrial, nSpike),0);
                Ca(nTrial,:) = Ca(nTrial,:) + exp(-Delta_t/tau_decay).*(1-exp(-Delta_t/tau_rise));
            end
        end
    end
    
    intNoiseCa = Ca + randn(size(Ca))*intNoise; 
    intNoiseCa(intNoiseCa<0) = 0;    
    CaImaging  = Fm./ (1+exp(-(Ca-Ca0)*beta));%Ca.^n ./ (K^n + Ca.^n);
%     intNoiseCaImaging = Fm* intNoiseCa.^n ./ (K^n + intNoiseCa.^n);
    intNoiseCaImaging = Fm./ (1+exp(-(intNoiseCa-Ca0)*beta));
    extNoiseCaImaging = CaImaging + randn(size(CaImaging))*extNoise; 

    
    if mod(rMean, 3) == 2
        subplot(1, 4, 2)
        hold on
        plot(timeSeriesData, z1(mean(CaImaging, 1)), ...
            '-', 'color', colorSet(rMean, :), 'linewid', 1);
        box off
        ylabel('normalized F')
        xlabel('Time (sec)')
        ylim([0 1])
        xlim([0 1])
        set(gca, 'TickDir', 'out')
        title ('no noise')

        subplot(1, 4, 3)
        hold on
        [~, idx] = max(mean(CaImaging, 1));
        min(idx)
        plot(timeSeriesData, z1(mean(intNoiseCaImaging, 1)), ...
            '-', 'color', colorSet(rMean, :), 'linewid', 1);
        box off
        ylabel('normalized F')
        xlabel('Time (sec)')
        ylim([0 1])  
        xlim([0 1])
        set(gca, 'TickDir', 'out')
        title ('internal noise')

        subplot(1, 4, 4)
        hold on
        plot(timeSeriesData, mean(extNoiseCaImaging, 1)/max(mean(CaImaging, 1)), ...
            '-', 'color', colorSet(rMean, :), 'linewid', 1);
        box off
        ylabel('normalized F')
        xlabel('Time (sec)')
        ylim([0 1])
        xlim([0 1])
        set(gca, 'TickDir', 'out')
        title ('external noise')
    end
    
end

setPrint(8*4, 6, 'S2CSuppRiseTime')
% close all