%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first p-value ramping time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffleConfounding.mat']);
% Gaussian filter for spiking data
sigma                         = 0.1 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 10;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);

if ~exist([PlotDir 'ConfoundingFactorFPVT'],'dir')
    mkdir([PlotDir 'ConfoundingFactorFPVT'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All area vs area of spiking recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nData          = 5;
load([TempDatDir DataSetList(nData).name '.mat'])
bumpStartPoint = getFPVT(nDataSet, filterInUse);
startTimePoint = DataSetList(nData).params.timeSeries(bumpStartPoint);
timeBins       = DataSetList(nData).params.timeSeries;
numT           = length(timeBins);
figure;
hold on;
pdfMat = ksdensity(startTimePoint(bumpStartPoint<numT), timeBins);
plot(timeBins, pdfMat, 'linewid', 2.0);

APLoc          = [DataSetList(nData).cellinfo.AP_axis];
MLLoc          = [DataSetList(nData).cellinfo.ML_axis];
pdfMat = ksdensity(startTimePoint(bumpStartPoint'<numT & APLoc>2400 & APLoc<2600 & MLLoc>1100 & MLLoc<1900), timeBins);
plot(timeBins, pdfMat, 'linewid', 2.0);
legend({'All', 'Sub'})
legend('location', 'south')
legend('boxoff')
gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 2.0)
hold off;
xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)])
xlabel('First P Val. time (ms)')
ylabel('Prob. Density')
set(gca, 'TickDir', 'out')
setPrint(8, 6, [PlotDir 'ConfoundingFactorFPVT/SingleUnitsFPVTAPML_' DataSetList(nData).name], 'svg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison across animals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmap = cbrewer('qual', 'Set1', 10, 'cubic');
DataSetToAnalysis = [1 6];
for nData          = DataSetToAnalysis
    load([TempDatDir DataSetList(nData).name '.mat'])
    bumpStartPoint = getFPVT(nDataSet, filterInUse);
    startTimePoint = DataSetList(nData).params.timeSeries(bumpStartPoint);
    timeBins       = DataSetList(nData).params.timeSeries;
    numT           = length(timeBins);
    [~, ~, anmIndex] = unique(cell2mat({DataSetList(nData).cellinfo.anmName}'), 'rows');
    figure;
    hold on;
    tColor = 0;
    
    for nAnm    = 1:anmIndex(end)
        if sum(anmIndex == nAnm) > 50
            tColor = tColor + 1;
            pdfMat = ksdensity(startTimePoint(bumpStartPoint<numT & anmIndex == nAnm), timeBins);
            plot(timeBins, pdfMat, 'linewid', 2.0, 'Color', cmap(tColor, :));
        end
    end
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)])
    xlabel('First P Val. time (ms)')
    ylabel('Prob. Density')
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'ConfoundingFactorFPVT/SingleUnitsFPVTANM_' DataSetList(nData).name], 'svg')
end

close all