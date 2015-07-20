%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Selectivity over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_3
%

addpath('../Func');
setDir;
load ([TempDatDir 'DataList.mat']);

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotZScoreImagesc(nDataSet, DataSetList(nData).params); 
    setPrint(8, 6, [PlotDir 'Single_Units_zScore_' DataSetList(nData).name], 'pdf')
end