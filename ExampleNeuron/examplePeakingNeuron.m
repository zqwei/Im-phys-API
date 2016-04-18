
function examplePeakingNeuron
    addpath('../Func');
    setDir;
    
    load ([TempDatDir 'DataListShuffle.mat']);
    nData = 1; % plot raster and psth
    load([TempDatDir DataSetList(nData).name '.mat'])
    spikeDataSet = nDataSet;   
    params = DataSetList(nData).params;
    nData = 4; % plot ca comparison
    load([TempDatDir DataSetList(nData).name '.mat'])
    caDataSet    = nDataSet;
    load ([TempDatDir 'DataListS2CModel.mat']);    
    nData = 1; % plot s2c model
    load([TempDatDir DataSetList(nData).name '.mat'])
    s2cDataSet   = nDataSet;
    
    meanCaDataSet = meanData(caDataSet);
    meanS2CDataSet = meanData(s2cDataSet);

%     peakingNeuronIndex = [97 244 378 441 459 236 354 387 457 492 239 370 436 ...
%         510 604 615 654 714 771 779 817 828 872 875 957 995];
    peakingNeuronIndex = 457;
    
    for nCell = peakingNeuronIndex
        figure;
        subplot(4, 1, 1)
        plotRaster(spikeDataSet(nCell), params);
        subplot(4, 1, 2)
        plotPSTH(spikeDataSet(nCell), params, 'Firing rate (Hz)');
        subplot(4, 1, 3)
        plotPSTH(s2cDataSet(nCell), params, 'DF/F');
        subplot(4, 1, 4)
        mCell = findSimilarCellToS2CModel(meanS2CDataSet(nCell,:), meanCaDataSet);
        plotPSTH(caDataSet(mCell), params, 'DF/F');
        setPrint(8, 6*4, [PlotDir 'SingleUnitsPeakFinderStats/SingleUnitsPeakExampleNeuron_' num2str(nCell, '%04d')])
        close all
    end

end

function plotRaster(spikeDataSet, params)
    spkTimes{1}    = spikeDataSet.unit_yes_trial_spk_time;
    spkTimes{2}    = spikeDataSet.unit_no_trial_spk_time;
    hold on;
    color_index    = [0.7  0 0; 0 0 0.7];
    spkLoc         = 0;
    for nPlot            = 1:2
        hold on;
        for nTrial       = 1:length(spkTimes{nPlot})
            spikeTimes   = spkTimes{nPlot}{nTrial};
            spikeTrial   = ones(length(spikeTimes), 1) * (nTrial + spkLoc);
            plot(spikeTimes, spikeTrial, '.', 'color', color_index(nPlot, :));
        end
        spkLoc     = spkLoc + length(spkTimes{nPlot}) + 3;
    end

    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0);
    hold off;
    ylim([1 spkLoc])
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    axis off
end

function plotPSTH(spikeDataSet, params, ylabelName)
    color_index    = [0.7  0 0; 0 0 0.7];
    hold on;
    nUnitData        = spikeDataSet.unit_yes_trial;
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
    nUnitData        = spikeDataSet.unit_no_trial;
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel(ylabelName);
    xlabel('Time (s)');
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    set(gca, 'TickDir', 'out')
end


function meanDataSet = meanData(nDataSet)
    meanDataSet = nan(length(nDataSet), size(nDataSet(1).unit_yes_trial, 2)*2);
    for nCell = 1:length(nDataSet)
        meanDataSet(nCell, :) = [mean(nDataSet(nCell).unit_yes_trial), mean(nDataSet(nCell).unit_no_trial)];
    end
end

function mIndex = findSimilarCellToS2CModel(nCellS2CDataSet, caDataSet)
    corrDat     = corr(nCellS2CDataSet', caDataSet');
    [~, mIndex] = max(corrDat);
    mIndex      = mIndex(1);
end