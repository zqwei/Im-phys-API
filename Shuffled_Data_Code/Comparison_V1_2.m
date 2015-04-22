%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Shape of activity distribution across trial periods (pre-sample,
%     sample, delay, response)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_2
% version 1.0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1 Histogram of the distribution of each units in different trial
%     periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% sampleSize            = 16;
% % m                     = ceil(sqrt(sampleSize));
% % numUnits              = 200;
% % sampleSeq             = randperm(numUnits);
% % sampleSeq             = sampleSeq(1:sampleSize);
% 
% load ('TempDat/DataList.mat');
% dist                  = {'Poisson', 'Normal', 'Normal', 'Normal'};
% 
% for nData             = 1:length(DataSetList)
%     load(['TempDat/' DataSetList(nData).name '.mat'])
%     plotHistActivityExampleUnits(nDataSet, sampleSeq, DataSetList(nData).params, dist{nData}); 
%     setPrint(m*4, m*3, ['Plot/Single_Units_Hist_' DataSetList(nData).name ], 'png')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.2 Histogram of the distribution of each units in different trial
%     periods (all neurons)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% load ('TempDat/DataList.mat');
% dist                  = {'Exp', 'Exp', 'Normal', 'Normal'};
% 
% for nData             = 1:length(DataSetList)
%     load(['TempDat/' DataSetList(nData).name '.mat'])
%     plotHistActivityPop(nDataSet, DataSetList(nData).params, dist{nData}); 
%     setPrint(m*4, m*3, ['Plot/Single_Units_Hist_' DataSetList(nData).name ], 'pdf')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_2
% version 2.0
%
%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
% dist                  = {'Poisson', 'Normal', 'Normal', 'Normal', 'Normal'};
% xlabels               = {'Firing rate(Hz)', 'dF/F', 'dF/F', 'dF/F', 'dF/F'};
% yAxes_set             = [-150 150; -100 100; -100 100; -200 200; -500 500];
% xAxes_set             = [0 20; -0.5 0.5; -0.4 0.41; -0.4 0.41; -0.4 0.41];
barSeries             = {0:1:60, -0.4:0.02:0.4, -0.4:0.02:0.4, -0.4:0.02:0.4, -0.4:0.02:0.4, -0.4:0.02:0.4};

% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotHistActivityPop(nDataSet, DataSetList(nData).params, dist{nData}, xlabels{nData},barSeries{nData},yAxes_set(nData,:),xAxes_set(nData,:)); 
%     m = 2;
%     setPrint(m*6, m*4.5, [PlotDir 'Single_Units_Hist_' DataSetList(nData).name ], 'pdf')
% end

nFactor = 1/2;
if ~exist([PlotDir '/Single_Units_Hist'],'dir')
    mkdir([PlotDir '/Single_Units_Hist'])
end


for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotHistActivityPopWithNonActiveNeurons(nDataSet(ephysCellIndex{nData}), DataSetList(nData).params, barSeries{nData}, nFactor); 
    m = 2;
    setPrint(m*8, m*6, [PlotDir 'Single_Units_Hist/Single_Units_Hist_' DataSetList(nData).name ], 'pdf')
end

close all