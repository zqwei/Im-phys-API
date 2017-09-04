%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single neuron choice probability index distribution (ROC), across
% trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

if ~exist([PlotDir 'SingleUnitsROC'],'dir')
    mkdir([PlotDir 'SingleUnitsROC'])
end

% 
% for nData             = 1:length(DataSetList)-1
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotROCPopAccLines(nDataSet, DataSetList(nData).params); 
%     set(gca, 'TickDir', 'out')
%     setPrint(8, 6, [PlotDir 'SingleUnitsROC/SingleUnitsROC_' DataSetList(nData).name])
% end
% 
% for nData             = 1:length(DataSetList)-1
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotROCPopAccLines(nDataSet(DataSetList(nData).ActiveNeuronIndex), DataSetList(nData).params); 
%     set(gca, 'TickDir', 'out')
%     setPrint(8, 6, [PlotDir 'SingleUnitsROC/SingleActUnitsROC_' DataSetList(nData).name])
% end
% 

% for nData             = 1:length(DataSetList)-1
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotROCPopLines(nDataSet(DataSetList(nData).ActiveNeuronIndex), DataSetList(nData).params); 
%     setPrint(8*3, 6, [PlotDir 'SingleUnitsROC/SingleActUnitsROCNoAcc_' DataSetList(nData).name])
% end

figure;
for nData             = [1 3 4]
    if nData == 1
        load([TempDatDir DataSetList(nData).name '.mat'])
        ActiveNeuronIndex = DataSetList(nData).ActiveNeuronIndex;
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
        ActiveNeuronIndex = DataSetList(nData).ActiveNeuronIndex(~neuronRemoveList);
%         load([TempDatDir DataSetList(nData).name '.mat'])
    end
    plotROCPopAccLinesSubplot(nDataSet(ActiveNeuronIndex), DataSetList(nData).params); 
    
end

setPrint(8*3, 6, [PlotDir 'SingleUnitsROC/SingleUnitsROC'])
% 
% close all