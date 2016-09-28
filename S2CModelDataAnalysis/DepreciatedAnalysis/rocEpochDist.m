%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single neuron choice probability index distribution (ROC), across
% trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CModel.mat']);

if ~exist([PlotDir 'S2CModel'],'dir')
    mkdir([PlotDir 'S2CModel'])
end



for nData             = [3 4]
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotROCPopLines(nDataSet(DataSetList(nData).ActiveNeuronIndex), DataSetList(nData).params); 
    setPrint(8*3, 6, [PlotDir 'S2CModel/SingleActUnitsROCNoAcc_' DataSetList(nData).name])
end


close all