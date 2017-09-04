%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single neuron choice probability index distribution (ROC), across
% trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rocEpochDistSummary
    addpath('../Func');
    setDir;
    load ([TempDatDir 'DataListS2CGP43Model.mat']);

    
    cmap = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

    for nData             = [2]
        if nData == 1
            load([TempDatDir 'Shuffle_Spikes.mat'])
%             ActiveNeuronIndex = DataSetList(nData).ActiveNeuronIndex;
        else
            load([TempDatDir DataSetList(nData).name '.mat'])
%             ActiveNeuronIndex = DataSetList(nData).ActiveNeuronIndex(~neuronRemoveList);
        end
%         nDataSet            = nDataSet(ActiveNeuronIndex);
        params              = DataSetList(nData).params;
        timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);     
        numPlots            = length(nDataSet);
        nBins               = 10;    
        histXout            = 0.5:0.02:1;
        histFreq            = zeros(length(histXout),length(timePoints) -1);

        for nPeriods        = 1:length(timePoints) -1
            nPeriodData     = dataInPeriods(nDataSet, timePoints, nPeriods);     
            areaInt         = zeros(numPlots, 1);
            for nPlot       = 1:numPlots
                nRocTData        = [nPeriodData(nPlot).unit_yes_trial; nPeriodData(nPlot).unit_no_trial];
                nRocOData        = [ones(size(nPeriodData(nPlot).unit_yes_trial)); zeros(size(nPeriodData(nPlot).unit_no_trial))];
                [tp, fp]         = RocFunc(nRocTData, nRocOData, nBins);
                areaInt(nPlot)   = intXY([tp, 1], [fp, 1]);
                areaInt(nPlot)   = max(areaInt(nPlot), 1-areaInt(nPlot));
            end
            histFreq(:,nPeriods) = ksdensity(areaInt, histXout);
            histFreq(:,nPeriods) = histFreq(:,nPeriods)/sum(histFreq(:,nPeriods));
    %         histFreq(:,nPeriods) = cumsum(histFreq(:,nPeriods))/sum(histFreq(:,nPeriods));
        end

        nTitles = {'pre-sample','sample','delay','response'};

        for nPeriods = 2:4
            subplot(1, 3, nPeriods-1)
            hold on
            plot(histXout, histFreq(:, nPeriods),'-', 'linewid', 1.0, 'color', cmap(nData, :));
            ylabel('% Accumulated Units')
            xlabel('Area under ROC');
            xlim([0.5 1])
            ylim([0 0.20])
            title(nTitles{nPeriods})
            box off
            hold off
        end
    end

    setPrint(8*3, 6, [PlotDir 'S2CGP43Model/SingleUnitsROC'])
    % 
    % close all
end

function areaInt       = intXY(vec_x, vec_y)
    
    areaInt            = sum((vec_y(1:end-1)+vec_y(2:end)).*(vec_x(2:end)-vec_x(1:end-1))/2);
    
end