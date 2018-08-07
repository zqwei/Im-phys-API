addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

if ~exist([PlotDir 'SingleUnitsExampleNeuronOutlier'],'dir')
    mkdir([PlotDir 'SingleUnitsExampleNeuronOutlier'])
end

color_index    = [0.7  0 0; 0 0 0.7];
thres          = 0.6;

for nData              = [3 4]
    load([TempDatDir DataSetList(nData).name '.mat'])
    params             = DataSetList(nData).params;
    for nUnit          = 1:length(nDataSet)
        nCellData      = nDataSet(nUnit);
        idx_yes_remove = mean(isnan(nCellData.unit_yes_trial_removeoutlier), 2)>thres;
        idx_no_remove  = mean(isnan(nCellData.unit_no_trial_removeoutlier), 2)>thres; 
        if sum(idx_yes_remove)>0 || sum(idx_no_remove)>0 
            figure;
            subplot(2, 2, 1)
            hold on;
            nUnitData        = nCellData.unit_yes_trial;
            stdYesTrial      = mean(var(nUnitData));
            stdYesTime       = var(mean(nUnitData));
            shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
                            std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
                            {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
            nUnitData        = nCellData.unit_no_trial;
            stdNoTrial       = mean(var(nUnitData));
            stdNoTime        = var(mean(nUnitData));
            shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
                std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
                {'k-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
            gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
            hold off;
            ylabel('DF/F');
            xlabel('Time (s)');
            xlim([params.timeSeries(1) params.timeSeries(end)]);
            set(gca, 'TickDir', 'out')
            title({'Before removal'; ...
                   ['ipsi: trial var:' num2str(stdYesTrial), ', temporal var:' num2str(stdYesTime)];
                   ['contra: trial var:' num2str(stdNoTrial), ', temporal var:' num2str(stdNoTime)]})

            subplot(2, 2, 3)
            hold on;
            nUnitData        = nCellData.unit_yes_trial(~idx_yes_remove, :);
            stdYesTrial      = mean(var(nUnitData));
            stdYesTime       = var(mean(nUnitData));
            shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
                            std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
                            {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
            nUnitData        = nCellData.unit_no_trial(~idx_no_remove, :);
            stdNoTrial       = mean(var(nUnitData));
            stdNoTime        = var(mean(nUnitData));
            shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
                std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
                {'k-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
            gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
            hold off;
            ylabel('DF/F');
            xlabel('Time (s)');
            xlim([params.timeSeries(1) params.timeSeries(end)]);
            set(gca, 'TickDir', 'out')
            
            title({'After removal'; ...
                   ['ipsi: trial var:' num2str(stdYesTrial), ', temporal var:' num2str(stdYesTime)];
                   ['contra: trial var:' num2str(stdNoTrial), ', temporal var:' num2str(stdNoTime)]})

            subplot(2, 2, 2)
            hold on;
            nUnitData        = nCellData.unit_yes_trial;
            if sum(idx_yes_remove)>0
                plot(params.timeSeries, nUnitData(idx_yes_remove, :), '-k','linewid', 1.0)
            end
            if sum(~idx_yes_remove)>0
                plot(params.timeSeries, nUnitData(~idx_yes_remove, :), '-','linewid', 1.0, 'color', color_index(1, :))
            end
            gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
            hold off;
            ylabel('DF/F');
            xlabel('Time (s)');
            xlim([params.timeSeries(1) params.timeSeries(end)]);
            set(gca, 'TickDir', 'out')
            title('Ipsi. trial')

            subplot(2, 2, 4)
            hold on;
            nUnitData        = nCellData.unit_no_trial;
            if sum(idx_no_remove)>0
                plot(params.timeSeries, nUnitData(idx_no_remove, :), '-k','linewid', 1.0)
            end
            if sum(~idx_no_remove)>0
                plot(params.timeSeries, nUnitData(~idx_no_remove, :), '-','linewid', 1.0, 'color', color_index(2, :))
            end
            gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
            hold off;
            ylabel('DF/F');
            xlabel('Time (s)');
            xlim([params.timeSeries(1) params.timeSeries(end)]);
            set(gca, 'TickDir', 'out')
            title('Contra. trial')
            
            setPrint(8*2, 6*2, [PlotDir 'SingleUnitsExampleNeuronOutlier/SingleUnitsExampleNeuronOutlier_' DataSetList(nData).name '_' num2str(nUnit, '%04d')])
            
            close all
        end  
    end
end