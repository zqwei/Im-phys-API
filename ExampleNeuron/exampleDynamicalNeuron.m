%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);


if ~exist([PlotDir 'SingleUnitsTscore'],'dir')
    mkdir([PlotDir 'SingleUnitsTscore'])
end

numExample             = 4;

color_index   = [0.7 0.0 0.0; % CL
                 0.0 0.0 0.7; % CR
                 1.0 0.7 0.7; % EL
                 0.7 0.7 1.0]; % ER

nData         = 1;
load([TempDatDir DataSetList(nData).name '.mat'])             
sigma                         = 0.15 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 11;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);  
textYLabels                   = {'Firing rate (Hz)', 'DF/F', 'DF/F', 'DF/F'};

for nData              = [1 3 4]
    load([TempDatDir 'LogPValueTscore_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
    load([TempDatDir DataSetList(nData).name '.mat'])             

    unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);

    figure;
    
    for nGroup         = 1:2
        neuronIndice   = find(unitGroup==nGroup);
        for nExample   = 1:numExample
            neuronIndex = neuronIndice(nExample);
            
            subplot(2, 4, nExample + (nGroup-1)*4)
            hold on;
            nUnitData        = nDataSet(neuronIndex).unit_yes_trial;
            smoothedUnitData = getGaussianPSTH (filterInUse, nUnitData, 2);
            shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
            nUnitData        = nDataSet(neuronIndex).unit_no_trial;
            smoothedUnitData = getGaussianPSTH (filterInUse, nUnitData, 2);
            shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
            gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
            hold off;
            ylabel(textYLabels{nData});
            xlabel('Time (s)');
            xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
        end
    end
    setPrint(8*4, 6*2, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreExampleNeuron_' DataSetList(nData).name], 'svg')
end