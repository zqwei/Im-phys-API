%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Single neuron choice probability index distribution (ROC), across
%     trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.1 Single neuron choice ROC curve, across trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_4

% sampleSize            = 64;
% m                     = ceil(sqrt(sampleSize));
% numUnits              = 100;
% sampleSeq             = randperm(numUnits);
% sampleSeq             = sampleSeq(1:sampleSize);
% 
% load ([TempDatDir 'DataListShuffle.mat']);
% % dist                  = {'Poisson', 'Normal', 'Normal', 'Normal'};
% 
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotROCExampleUnits(nDataSet, sampleSeq, DataSetList(nData).params); 
%     setPrint(m*4, m*3, [PlotDir 'Single_Units_ROC_' DataSetList(nData).name], 'pdf')
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.2 ROC curve, across trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_4


% load ('TempDat/DataList.mat');
% 
% for nData             = 1:length(DataSetList)
%     load(['TempDat/' DataSetList(nData).name '.mat'])
%     plotROCPop(nDataSet, DataSetList(nData).params); 
%     setPrint(2*4, 2*3, ['Plot/Single_Units_ROC_' DataSetList(nData).name], 'pdf')
% end


% addpath('../Func');
% setDir;
% load ([TempDatDir 'DataListShuffle.mat']);
% 
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotROCPop(nDataSet, DataSetList(nData).params); 
%     setPrint(2*4, 2*3, [PlotDir 'Single_Units_ROC_' DataSetList(nData).name], 'pdf')
% end

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

if ~exist([PlotDir '/Single_Units_ROC'],'dir')
    mkdir([PlotDir '/Single_Units_ROC'])
end

% for nData             = fileToAnalysis
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotROCPopAccLines(nDataSet, DataSetList(nData).params); 
%     setPrint(8, 6, [PlotDir 'Single_Units_ROC/Single_Units_ROC_' DataSetList(nData).name], 'pdf')
% end

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotROCPopAccLines(nDataSet, DataSetList(nData).params); 
%     title(strrep(DataSetList(nData).name, '_', ' '))
    setPrint(8, 6, [PlotDir 'Single_Units_ROC/SingleUnitsROC_' DataSetList(nData).name], 'pdf')
end

for nData             = 2:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotROCPopAccLines(nDataSet(DataSetList(nData).ActiveNeuronIndex), DataSetList(nData).params); 
%     title(strrep(DataSetList(nData).name, '_', ' '))
    setPrint(8, 6, [PlotDir 'Single_Units_ROC/SingleActUnitsROC_' DataSetList(nData).name], 'pdf')
end


for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotROCPopLines(nDataSet(DataSetList(nData).ActiveNeuronIndex), DataSetList(nData).params); 
%     title(strrep(DataSetList(nData).name, '_', ' '))
    setPrint(8, 6, [PlotDir 'Single_Units_ROC/SingleActUnitsROCNoAcc_' DataSetList(nData).name], 'pdf')
end


for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotPVALUEPopLines(nDataSet(DataSetList(nData).ActiveNeuronIndex), DataSetList(nData).params); 
%     title(strrep(DataSetList(nData).name, '_', ' '))
    setPrint(8, 6, [PlotDir 'Single_Units_ROC/SingleActUnitsPVALUENoAcc_' DataSetList(nData).name], 'pdf')
end

close all


% for nData             = fileToAnalysis
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotROCPopCompLines(nDataSet, DataSetList(nData).params); 
%     setPrint(8, 6, [PlotDir 'Single_Units_ROC/Single_Units_PreSampleVSampleROC_' DataSetList(nData).name], 'pdf')
% end

close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 4.2 % neuron with strong ROC across trial
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Comparison_V1_4
% 
% 
% addpath('../Func');
% setDir;
% load ([TempDatDir 'DataListShuffle.mat']);
% nColor = {'r', 'b', 'g', 'y', 'k'};
% ROCthres = 0.6;
% figure;
% 
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotROCPopTime(nDataSet, DataSetList(nData).params, ROCthres, nColor{nData}); 
% end
% 
% setPrint(2*4, 2*3, [PlotDir 'Single_Units_ROC_' DataSetList(nData).name], 'pdf')