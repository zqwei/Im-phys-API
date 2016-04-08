%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single neuron choice probability index distribution (ROC), across
% trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);

if ~exist([PlotDir 'ModeledSingleUnitsROC'],'dir')
    mkdir([PlotDir 'ModeledSingleUnitsROC'])
end


for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotROCPopAccLines(nDataSet, DataSetList(nData).params); 
    setPrint(8, 6, [PlotDir 'ModeledSingleUnitsROC/SingleUnitsROC_' DataSetList(nData).name], 'pdf')
end

% for nData             = 1:length(DataSetList)-1
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotROCPopAccLines(nDataSet(DataSetList(nData).ActiveNeuronIndex), DataSetList(nData).params); 
%     setPrint(8, 6, [PlotDir 'ModeledSingleUnitsROC/SingleActUnitsROC_' DataSetList(nData).name], 'pdf')
% end
% 
% 
% for nData             = 1:length(DataSetList)-1
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotROCPopLines(nDataSet(DataSetList(nData).ActiveNeuronIndex), DataSetList(nData).params); 
%     setPrint(8, 6, [PlotDir 'ModeledSingleUnitsROC/SingleActUnitsROCNoAcc_' DataSetList(nData).name], 'pdf')
% end


close all