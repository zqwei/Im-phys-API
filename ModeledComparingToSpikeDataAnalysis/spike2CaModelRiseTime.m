%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% low firing activity to calcium analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Fm              = 21.3240;
K               = 13.9248;
n               = 1.7531;
tau_rise        = 0.0728;
tau_decay       = 1.4551;

timeSeriesData  = 0:0.01:5;
numTrial        = 20;
intNoise        = 2.5;
extNoise        = 2.5;

colorSet        = cool(20);

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
    

    subplot(3, 3, 1)
    hold on;
    for nTrial = 1:numTrial
        plot(preSamplePoints(nTrial, :), ones(1, 100)*nTrial + rMean*(numTrial+5), '.', 'color', colorSet(rMean, :))
    end

    text(0, rMean*(numTrial+5)+13, [num2str(rMean) ' Hz'])

    xlim([timeSeriesData(1) timeSeriesData(end)])
    ylim([numTrial+5 rMean*(numTrial+5)+numTrial+5+0.5])
    axis off

    subplot(2, 3, 4)
    hold on
    plot(timeSeriesData, mean(CaImaging, 1), '-', 'color', colorSet(rMean, :), 'linewid', 2);
    box off
    ylabel('F')
    xlabel('Time (sec)')
    
    subplot(2, 3, 3)
    hold on
    if ~isempty(find(mean(CaImaging, 1)>frThres, 1))
        plot(rMean, timeSeriesData(find(mean(CaImaging, 1)>frThres, 1)), 'bo', 'linewid', 2)
    end
    box off
    ylabel('Rising time to F=10 (sec)')
    xlabel('Spiking rate (Hz)')
    
    
    subplot(2, 3, 5)
    hold on
    plot(timeSeriesData, mean(intNoiseCaImaging, 1), '-', 'color', colorSet(rMean, :), 'linewid', 2);
    box off
    ylabel('F')
    xlabel('Time (sec)')
    
    subplot(2, 3, 3)
    hold on
    if ~isempty(find(mean(CaImaging, 1)>frThres, 1))
        plot(rMean, timeSeriesData(find(mean(CaImaging, 1)>frThres, 1)), 'ko', 'linewid', 2)
    end
    box off
    ylabel('Rising time to F=10 (sec)')
    xlabel('Spiking rate (Hz)')
    
    
    subplot(2, 3, 6)
    hold on
    plot(timeSeriesData, mean(extNoiseCaImaging, 1), '-', 'color', colorSet(rMean, :), 'linewid', 2);
    box off
    ylabel('F')
    xlabel('Time (sec)')
    
    subplot(2, 3, 3)
    hold on
    if ~isempty(find(mean(CaImaging, 1)>frThres, 1))
        plot(rMean, timeSeriesData(find(mean(CaImaging, 1)>frThres, 1)), 'ro', 'linewid', 2)
    end
    box off
    ylabel(['Rising time to F=' num2str(frThres) ' (sec)'])
    xlabel('Spiking rate (Hz)')
    
    
    
end

% setPrint(8*3, 8, 'Spike2CaAnalysis', 'pdf')