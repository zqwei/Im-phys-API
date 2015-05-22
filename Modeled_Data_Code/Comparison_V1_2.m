%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Shape of activity distribution across trial periods (pre-sample,
%     sample, delay, response)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_2
% version 1.0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_2
% version 2.0
%
%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);
fileToAnalysis                = [1 4 5]; % 1:length(DataSetList);

% dist                  = {'Poisson', 'Normal', 'Normal', 'Normal', 'Normal'};
xlabels               = {'dR/R', 'dF/F', 'dF/F', 'dF/F', 'dF/F', 'dF/F'};
% yAxes_set             = [-150 150; -100 100; -100 100; -200 200; -500 500];
% xAxes_set             = [0 20; -0.5 0.5; -0.4 0.41; -0.4 0.41; -0.4 0.41];
barSeries             = {-2:0.04:2, -2:0.04:2, -2:0.04:2, -2:0.04:2, -2:0.04:2, -2:0.04:2};

% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotHistActivityPop(nDataSet, DataSetList(nData).params, dist{nData}, xlabels{nData},barSeries{nData},yAxes_set(nData,:),xAxes_set(nData,:)); 
%     m = 2;
%     setPrint(m*6, m*4.5, [PlotDir 'Single_Units_Hist_' DataSetList(nData).name ], 'pdf')
% end

nFactor = 1/7;
if ~exist([PlotDir '/SingleModel_Units_Hist'],'dir')
    mkdir([PlotDir '/SingleModel_Units_Hist'])
end


for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotHistActivityPopWithNonActiveNeurons(nDataSet, DataSetList(nData).params, barSeries{nData}, nFactor, xlabels{nData}); %[DataSetList(nData).cellinfo(:).cellType]==1
    m = 2;
    setPrint(m*8, m*6, [PlotDir 'SingleModel_Units_Hist/Single_Units_Hist_' DataSetList(nData).name ], 'pdf')
end


close all