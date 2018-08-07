function exampleS2CNeuron
    addpath('../Func');
    setDir;
    
    load ([TempDatDir 'DataListShuffle.mat']);
    nData = 1; % plot raster and psth
    load([TempDatDir DataSetList(nData).name '_old.mat'])
    spikeDataSet = nDataSet;   
    
    load('ParamsFitCells_S2CModel_Sim.mat');
    S2Cparams   = params;
    clear params;
    
    params.frameRate       =  29.68/2;
    params.binsize         =  1/params.frameRate;
    params.polein          =  -2.6;
    params.poleout         =  -1.3;
    minTimeToAnalysis      =  round(-3.1 * params.frameRate);
    maxTimeToAnalysis      =  round(2.0 * params.frameRate);
    params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
    params.timeSeries      = params.timeWindowIndexRange * params.binsize;
    params.expression      = 'None';
    idTauD      = 1;
    idN         = 1;
    idK         = 1;

    % spikeDataSet  = nDataSet([491, 189, 401, 555, 55]);
    spikeDataSet  = nDataSet([219, 568, 970]);
    params.Fm              = 21.3240*ones(size(spikeDataSet));
    params.K               = S2Cparams(idN).K*ones(size(spikeDataSet));
    params.n               = S2Cparams(idN).n*ones(size(spikeDataSet));
    params.tau_r           = S2Cparams(1).tau_r*ones(size(spikeDataSet));
    params.tau_d           = S2Cparams(idTauD).tau_d*ones(size(spikeDataSet));%1.1*ones(size(spikeDataSet));%S2Cparams(idTauD).tau_d*ones(size(spikeDataSet));
    params.intNoise        = 0.5;
    params.extNoise        = 1.5;
    load('FineTunedNLParams.mat');
    g = @(p,x) p(1) + p(2)./ (1 + 10.^((p(3)-x)*p(4)));
    ext_noise = 0.15;
    numTrial = 50;
    % nlParams   = nlParams([491, 189, 401, 555, 55], :);
    nlParams   = nlParams([219, 568, 970], :);
    
    slopepara  = [0.1 0.3 0.5 0.8 1.5 2];
    tau_d_set = [0.28, 0.5, 1.1, 1.7, 2.5, 3.1];
    
    for nDecay  = 1:6
        params.tau_d           = 1.7*ones(size(spikeDataSet));%tau_d_set(nDecay)*ones(size(spikeDataSet));
        nDataSet               = getFakeCaImagingData(spikeDataSet, params);
        for nUnit  = 1:length(nDataSet)
            param  = nlParams(nUnit, :);
            param(1) = 0;
            param(2) = 3;
            param(3) = 1/slopepara(nDecay) * param(3);
            param(4) = slopepara(nDecay) * param(4);
            yesNoise = randn(numTrial, 77)*ext_noise; %squeeze(noiseRates{nData}(nUnit, 1, :, :))';
            noNoise  = randn(numTrial, 77)*ext_noise; %squeeze(noiseRates{nData}(nUnit, 2, :, :))';
            nDataSet(nUnit).unit_yes_trial = g(param, nDataSet(nUnit).unit_yes_trial_linear+...
                                                params.intNoise*randn(size(nDataSet(nUnit).unit_yes_trial_linear))) + ...
                                                yesNoise(randpermLargeK(numTrial, size(nDataSet(nUnit).unit_yes_trial, 1)), :);
            nDataSet(nUnit).unit_no_trial  = g(param, nDataSet(nUnit).unit_no_trial_linear+...
                                                params.intNoise*randn(size(nDataSet(nUnit).unit_no_trial_linear))) + ...
                                                noNoise(randpermLargeK(numTrial, size(nDataSet(nUnit).unit_no_trial, 1)), :);
            nDataSet(nUnit).unit_yes_error = g(param, nDataSet(nUnit).unit_yes_error_linear) + ...
                                                yesNoise(randpermLargeK(numTrial, size(nDataSet(nUnit).unit_yes_error, 1)), :);
            nDataSet(nUnit).unit_no_error  = g(param, nDataSet(nUnit).unit_no_error_linear) + ...
                                                noNoise(randpermLargeK(numTrial, size(nDataSet(nUnit).unit_no_error, 1)), :);
        end



%         params = DataSetList(nData).params;

        for nCell = 1:length(nDataSet)
            subplot(length(nDataSet), 6, (nCell-1)*6+nDecay)
            plotPSTH(nDataSet(nCell), params, 'dF/F')
        end
    end

    setPrint(8*6, 6*4, [PlotDir 'SingleUnitsTscore/SingleUnitsS2CExampleNeuron'])
%     close all
end


function plotPSTH(spikeDataSet, params, ylabelName)
    sigma                         = 0.15 / params.binsize; % 200 ms
    filterLength                  = 11;
    filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
    filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
    filterInUse                   = filterInUse / sum (filterInUse); 

    color_index    = [0 0 0.7; 0.7  0 0];
    hold on;
    nUnitData        = spikeDataSet.unit_yes_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    meanBaseline     = mean([spikeDataSet.unit_yes_trial; spikeDataSet.unit_no_trial]);
    meanBaseline     = mean(meanBaseline(1:8));
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1)-meanBaseline,...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
    nUnitData        = spikeDataSet.unit_no_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1)-meanBaseline,...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel(ylabelName);
    xlabel('Time (s)');
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    set(gca, 'TickDir', 'out')
end