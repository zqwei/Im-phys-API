%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Single neuron choice probability index distribution (ROC), across
%     trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_4



addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);

if ~exist([PlotDir '/SingleModel_Units_ROC'],'dir')
    mkdir([PlotDir '/SingleModel_Units_ROC'])
end

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotROCPopAccLines(nDataSet, DataSetList(nData).params); 
    setPrint(8, 6, [PlotDir 'SingleModel_Units_ROC/Single_Units_ROC_' DataSetList(nData).name], 'pdf')
end

close all