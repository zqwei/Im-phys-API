
function examplePeakingNeuron
    addpath('../Func');
    setDir;
    
    load ([TempDatDir 'DataListShuffle.mat']);
    nData = 1; % plot raster and psth
    load([TempDatDir DataSetList(nData).name '_old.mat'])
    spikeDataSet = nDataSet;   
    params = DataSetList(nData).params;
    dynamicalNeuronIndex = [1047 982 972 810 751 742 728 555 401 395 307 264 221 161 137 1041];
    dynamicalNeuronIndex = dynamicalNeuronIndex([6, 5, 4, 8, 16 ]);
    params.Fm = 1;
    params.K  = 1;
    params.n  = 1;
    params.tau_r = 0.060;
    params.tau_d = 1.7;
    params.intNoise = 0;
    params.extNoise = 0;
    ext_noise = params.extNoise;
    g = @(p,x) p(1) + p(2)./ (1 + 10.^((p(3)-x)*p(4)));
    param = [0, 1, 0.2, 3];
    
    numTrial    = 100;
    
    for nCellid = 1:length(dynamicalNeuronIndex)
        nCell   = dynamicalNeuronIndex(nCellid);
        nCellDataSet  = getFakeCaImagingData(spikeDataSet(nCell), params);
        yesNoise = randn(numTrial, 77)*ext_noise; %squeeze(noiseRates{nData}(nUnit, 1, :, :))';
        noNoise  = randn(numTrial, 77)*ext_noise; %squeeze(noiseRates{nData}(nUnit, 2, :, :))';
        mean_pre_sample             = max(mean(nCellDataSet.unit_yes_trial_linear(:,1:8)));
        param(3)                    = mean_pre_sample;
        nCellDataSet.unit_yes_trial = nCellDataSet.unit_yes_trial_linear-mean_pre_sample;
        nCellDataSet.unit_yes_trial = nCellDataSet.unit_yes_trial/max(mean(nCellDataSet.unit_yes_trial));
        nCellDataSet.unit_no_trial  = g(param, nCellDataSet.unit_no_trial_linear+...
                                            params.intNoise*randn(size(nCellDataSet.unit_no_trial_linear))) + ...
                                            noNoise(randpermLargeK(numTrial, size(nCellDataSet.unit_no_trial, 1)), :);
        nCellDataSet.unit_yes_error = g(param, nCellDataSet.unit_yes_error_linear) + ...
                                            yesNoise(randpermLargeK(numTrial, size(nCellDataSet.unit_yes_error, 1)), :);
        nCellDataSet.unit_no_error  = g(param, nCellDataSet.unit_no_error_linear) + ...
                                            noNoise(randpermLargeK(numTrial, size(nCellDataSet.unit_no_error, 1)), :);
        
        subplot(length(dynamicalNeuronIndex), 2, 1+nCellid*2-2)
        plotPSTH(spikeDataSet(nCell), params, 'Firing rate (Hz)');

        subplot(length(dynamicalNeuronIndex), 2, 2+nCellid*2-2)
        plotPSTH(nCellDataSet, params, 'DF/F');        
    end
    setPrint(8*2, 6*length(dynamicalNeuronIndex'), [PlotDir 'SingleUnitsImagescWithSort/SingleUnitsTscoreExampleNeuronS2C_Linear'])
%     close all

end

function plotRaster(spikeDataSet, params)
    spkTimes{1}    = spikeDataSet.unit_yes_trial_spk_time;
    spkTimes{2}    = spikeDataSet.unit_no_trial_spk_time;
    hold on;
    color_index    = [0 0 0.7; 0.7  0 0];
    spkLoc         = 0;
    for nPlot            = 1:2
        hold on;
        for nTrial       = 1:20
            spkLoc       = spkLoc + 1;
            spikeTimes   = spkTimes{nPlot}{nTrial};
            spikeTrial   = ones(length(spikeTimes), 1) * spkLoc;
            plot(spikeTimes, spikeTrial, '.', 'color', color_index(nPlot, :));
        end
        spkLoc     = spkLoc + 1;%length(spkTimes{nPlot}) + 3;
    end

    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0);
    hold off;
    ylim([1 41])
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    axis off
end

function plotDff(spikeDataSet, params)
%     sigma                         = 0.15 / params.binsize; % 200 ms
%     filterLength                  = 11;
%     filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
%     filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
%     filterInUse                   = filterInUse / sum (filterInUse); 

    nUnitData        = spikeDataSet.unit_yes_trial;
%     yesUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
    yesUnitData      = nUnitData;
    nUnitData        = spikeDataSet.unit_no_trial;
    noUnitData       = nUnitData;
%     noUnitData       = getGaussianPSTH (filterInUse, nUnitData, 2);
    cmin             = min([mean(yesUnitData), mean(noUnitData)]);
    cmax             = max([mean(yesUnitData), mean(noUnitData)]);
    


    hold on
    actMat = nan(41, size(spikeDataSet.unit_yes_trial, 2));
    actMat(1:20, :) = yesUnitData(1:20, :);
    actMat(22:end, :) = noUnitData(1:20, :);
    imagesc(params.timeSeries, 1:41, actMat);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','w','Linestyle','--','linewid', 1.0);
    hold off;
    ylim([1 41])
    colormap(gray)
    
    caxis([cmin, cmax]);
%     colorbar
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    axis off
end


function plotPSTH(spikeDataSet, params, ylabelName)
    sigma                         = 0.10 / params.binsize; % 200 ms
    filterLength                  = 11;
    filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
    filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
    filterInUse                   = filterInUse / sum (filterInUse); 

    color_index    = [0 0 0.7; 0.7  0 0];
    hold on;
    nUnitData        = spikeDataSet.unit_yes_trial;    
%     nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
    x = mean(nUnitData, 1);
    [max_v, maxid]  = max(mean(nUnitData, 1));
    plot(params.timeSeries(maxid), max_v, 'go');
    nUnitData        = spikeDataSet.unit_yes_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    plot(params.timeSeries, mean(nUnitData, 1), '-k');
    maxid
    y = mean(nUnitData, 1);
    x(maxid)/y(maxid)
%     shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
%         std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
%         {'k-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
%     ylabel(ylabelName);
%     xlabel('Time (s)');
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