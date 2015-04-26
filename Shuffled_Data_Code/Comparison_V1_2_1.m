%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Shape of activity distribution across trial periods (pre-sample,
%     sample, delay, response)
%     negative neurons in Figure 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_2
% version 1.0

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffleRaw.mat']);
DataSetListRaw       = DataSetList;
load ([TempDatDir 'DataListShuffle.mat']);
fileToAnalysis                = [4 5]; % 1:length(DataSetList);

% dist                  = {'Poisson', 'Normal', 'Normal', 'Normal', 'Normal'};
xlabels               = {'Firing rate(Hz)', 'dF/F', 'dF/F', 'dF/F', 'dF/F'};
barSeries             = {0:1:60, -0.6:0.02:0.6, -0.6:0.02:0.6, -0.6:0.02:0.6, -0.6:0.02:0.6, -0.6:0.02:0.6};

nFactor = 1/2;
if ~exist([PlotDir '/Single_Units_Hist'],'dir')
    mkdir([PlotDir '/Single_Units_Hist'])
end


for nData             = fileToAnalysis
    load([TempDatDir DataSetListRaw(nData).name '.mat'])
    nDataSetRaw       = nDataSet;
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotNegativeActiveNeurons(nDataSet(ephysCellIndex{nData}), nDataSetRaw(ephysCellIndex{nData}), DataSetList(nData).params, xlabels{nData}); 
    setPrint(2*8, 9*6, [PlotDir 'Single_Units_Hist/Single_Units_NegativeNeurons_' DataSetList(nData).name ], 'pdf')
end

close all