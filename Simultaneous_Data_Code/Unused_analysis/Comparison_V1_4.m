%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Single neuron choice probability index distribution (ROC), across
%     trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_4

addpath('../Func');
setDir;
load ([TempDatDir 'DataList.mat']);

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotROCPop(nDataSet, DataSetList(nData).params); 
    setPrint(2*4, 2*3, [PlotDir 'Single_Units_ROC__' DataSetList(nData).name], 'pdf')
end
