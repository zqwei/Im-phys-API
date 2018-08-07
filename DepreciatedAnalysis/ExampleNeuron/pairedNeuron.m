neuronIndex          = [504 219];
figure;
for nCell            = 1:length(neuronIndex)
    nCellid          = neuronIndex(nCell);
    subplot(length(neuronIndex), 2, 2*(nCell-1)+1)
    hold on;
    nUnitData        = spikeDataSet(nCellid).unit_yes_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
    nUnitData        = spikeDataSet(nCellid).unit_no_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('Spikes /s');
    xlabel('Time (s)');
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    set(gca, 'TickDir', 'out')
    subplot(length(neuronIndex), 2, 2*(nCell-1)+2)
    hold on;
    nUnitData        = s2cDataSet(nCellid).unit_yes_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
    nUnitData        = s2cDataSet(nCellid).unit_no_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
    ylim([0 15])
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('df/f');
    xlabel('Time (s)');
    
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    set(gca, 'TickDir', 'out')
end

setPrint(8*2, 6*2, 'Similar')