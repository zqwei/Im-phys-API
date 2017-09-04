addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

if ~exist([PlotDir 'SingleUnitsExampleNeuronAfterOutlierRemoval'],'dir')
    mkdir([PlotDir 'SingleUnitsExampleNeuronAfterOutlierRemoval'])
end

color_index    = [0.7  0 0; 0 0 0.7];
thres          = 0.6;

for nData              = [3 4]
    load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    params             = DataSetList(nData).params;
    for nUnit          = 1:50;%1:length(nDataSet)
        nCellData      = nDataSet(nUnit);
        figure;
        hold on;
        nUnitData        = nCellData.unit_yes_trial;
        plot(params.timeSeries, nCellData.unit_yes_trial, '-','linewid', 1.0, 'color', color_index(1, :))
        plot(params.timeSeries, nCellData.unit_no_trial, '-','linewid', 1.0, 'color', color_index(2, :))
        gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
        hold off;
        ylabel('DF/F');
        xlabel('Time (s)');
        xlim([params.timeSeries(1) params.timeSeries(end)]);
        set(gca, 'TickDir', 'out')
        title('Contra. trial')
        setPrint(8*2, 6, [PlotDir 'SingleUnitsExampleNeuronAfterOutlierRemoval/SingleUnitsExampleNeuronAfterOutlierRemoval' DataSetList(nData).name '_' num2str(nUnit, '%04d')])
        close all
    end
end