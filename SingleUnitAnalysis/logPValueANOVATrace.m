%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% AP cell                : sensory
% LR cell                : decision
% CE cell                : reward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);


if ~exist([PlotDir 'SingleUnitsAnova'],'dir')
    mkdir([PlotDir 'SingleUnitsAnova'])
end

cmap = cbrewer('qual', 'Set1', 10, 'cubic');

thresLogP     = -log(0.01);
ylimSet       = [0.2, 0, 0.03, 0.1];

for nData      = [1 3 4]
    if nData == 1
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end
    sigma                         = 0.10 / DataSetList(1).params.binsize; % 200 ms
    filterLength                  = 10;
    filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
    filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
    filterInUse                   = filterInUse / sum (filterInUse);
    logPValue     = getLogPValueSpike(nDataSet, filterInUse);
    
    
    percentageNeuron = nan(size(logPValue, 3), 4);
    
    for nTime     = 1:size(logPValue, 3)
        unitGroup = squeeze(logPValue(:, :, nTime));
        gAP       = unitGroup(:, 1) > thresLogP;
        gLR       = unitGroup(:, 2) > thresLogP;
        gCE       = unitGroup(:, 3) > thresLogP;
        
        sumSelective = sum(gAP | gLR | gCE);
        sumUnit      = length(gAP);
        sumAP        = sum(gAP & ~gLR & ~gCE);
        sumLR        = sum(~gAP &  gLR & ~gCE);
        sumCE        = sum(~gAP & ~gLR &  gCE);
        percentageNeuron(nTime, 1) = sumAP/sumUnit;
        percentageNeuron(nTime, 2) = sumLR/sumUnit;
        percentageNeuron(nTime, 3) = sumCE/sumUnit;
        percentageNeuron(nTime, 4) = (sumSelective - (sumAP+sumLR+sumCE))/sumUnit;
    end
    
    figure;
    for nType  = 1:4
        subplot(4, 1, nType)
        hold on
        plot(DataSetList(nData).params.timeSeries, percentageNeuron(:, nType), '-k','linewid', 2.0)
        ylim([0 ylimSet(nData)])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
        xlim([params.timeSeries(1) params.timeSeries(end)]);
        xlabel('Time (s)');
        ylabel('Frac. cell')
        hold off
        box off
    end
    
    setPrint(8, 6*4, [PlotDir 'SingleUnitsAnova/SingleUnitsAnovaTime_' DataSetList(nData).name])
end


