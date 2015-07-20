%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Shape of activity distribution across trial periods (pre-sample,
%     sample, delay, response)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_2
%
addpath('../Func');
setDir;
load ([TempDatDir 'DataList.mat']);
dist                  = {'Poisson', 'Poisson', 'Normal', 'Normal', 'Normal'};
xlabels               = {'Firing rate(Hz)', 'dF/F', 'dF/F', 'dF/F', 'dF/F'};

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotHistActivityPop(nDataSet, DataSetList(nData).params, dist{nData}, xlabels{nData}); 
    m = 2;
    setPrint(m*6, m*4.5, [PlotDir 'Single_Units_Hist_' DataSetList(nData).name ], 'pdf')
end