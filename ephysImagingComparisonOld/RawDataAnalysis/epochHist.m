% Examine the firing rates in ephys vs whole cell recording data

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

% ephys
nData    = 10;
params   = DataSetList(nData).params;
timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
load([TempDatDir DataSetList(nData).name '.mat'])
m        = 2;

figure;
barSeries   = 0:0.3:15;
dist        = 'Poisson';
for nPlot   = 1:4
    nPeriodData = dataInPeriods(nDataSet, timePoints, nPlot);
    subplot(m, m, nPlot)
    hold on;
    barData     = nan(length(nPeriodData), 1);
    for nUnit   = 1:length(nPeriodData)
        barData(nUnit) = mean(nPeriodData(nUnit).unit_yes_trial);
    end
    barSign     = 1;
    barHistWithDist(barData(:), dist, '', barSeries, 'b', barSign); 
    barData     = nan(length(nPeriodData), 1);
    for nUnit   = 1:length(nPeriodData)
        barData(nUnit) = mean(nPeriodData(nUnit).unit_no_trial);
    end
    barSign     = -1;
    barHistWithDist(barData(:), dist, '', barSeries, 'r', barSign); 
    hold off;
end

% % whole cell
% nData = 9;
% load([TempDatDir DataSetList(nData).name '.mat'])
% for nPlot   = 1:4
%     nPeriodData = dataInPeriods(nDataSet, timePoints, nPlot);
%     subplot(m, m, nPlot)
%     hold on;
%     barData     = nan(length(nPeriodData), 1);
%     for nUnit   = 1:length(nPeriodData)
%         barData(nUnit) = mean(nPeriodData(nUnit).unit_yes_trial);
%     end
%     barSign     = 1;
%     barHistWithDist(barData(:), dist, '', barSeries, 'k', barSign); 
%     barData     = nan(length(nPeriodData), 1);
%     for nUnit   = 1:length(nPeriodData)
%         barData(nUnit) = mean(nPeriodData(nUnit).unit_no_trial);
%     end
%     barSign     = -1;
%     barHistWithDist(barData(:), dist, '', barSeries, 'g', barSign); 
%     hold off;
% end