%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12. GPFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_12

addpath('../Func');
setDir;
load ([TempDatDir 'DataList.mat']);
addNoise         = [1 1 0 0 0];
numTrials           = 400;
numTestTrials       = 100;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = rand(numTrainingTrials, 1) > 0.5;
testTargets         = rand(numTestTrials, 1) > 0.5;
totTargets          = [testTargets; trainingTargets];

addpath('../../../../../Empirical_Data_Analysis_Code/GPFA')
addpath('../../../../../Empirical_Data_Analysis_Code/GPFA/core_gpfa')
addpath('../../../../../Empirical_Data_Analysis_Code/GPFA/core_twostage')
addpath('../../../../../Empirical_Data_Analysis_Code/GPFA/plotting')
addpath('../../../../../Empirical_Data_Analysis_Code/GPFA/util')
addpath('../../../../../Empirical_Data_Analysis_Code/GPFA/util/precomp')
addpath('../../../../../Empirical_Data_Analysis_Code/GPFA/util/invToeplitz')

minNumUnit          = 7;
numFolds            = 4;
runIdx              = 0;
kernSD              = 30;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12.1 GPFA for single session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for nData            = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])  
%     numUnit          = arrayfun(@(x) length(x.nUnit), nDataSet3D, 'UniformOutput', false);
%     numUnit          = [numUnit{:}];
%     [numUnit, sortUnit]  = sort(numUnit,'descend');    
%     nDataSet3D       = nDataSet3D(sortUnit);
%     numSession       = sum(numUnit>=minNumUnit);
%     for nPlot        = 1:numSession
%         dat          = getGPFAData3D(nDataSet3D(nPlot));
%         runIdx       = runIdx + 1;
%         xDim         = ceil(numUnit(nPlot)/3*2);
%         neuralTraj(runIdx, dat, 'method', 'gpfa', 'xDim', xDim, 'numFolds', numFolds);
%         method       = plotPredErrorVsDim(runIdx, kernSD,'plotOn',false);
%         [~, xDim]    = min(method(5).sse);
%         result       = neuralTraj(runIdx, dat, 'method', 'gpfa', 'xDim', xDim);
%         [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
%         plotEachDimVsTimeGPFA(seqTrain, DataSetList(nData).params);
%         setPrint(4*8, ceil(xDim/4)*6, [PlotDir 'GPFA_' DataSetList(nData).name '_Session_' num2str(nPlot)], 'tif')
%     end
% end
% 
% close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12.2 GPFA for collected data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% runIdx               = 14;
% 
% for nData            = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     dat              = getGPFAData(nDataSet, totTargets);
%     runIdx           = runIdx + 1;
%     xDim             = ceil(size(dat(1).spikes,1)/3*2);
%     neuralTraj(runIdx, dat, 'method', 'gpfa', 'xDim', xDim, 'numFolds', numFolds);
%     method           = plotPredErrorVsDim(runIdx, kernSD,'plotOn',false);
%     [~, xDim]        = min(method(5).sse);
%     result           = neuralTraj(runIdx, dat, 'method', 'gpfa', 'xDim', xDim);
%     [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
%     plotEachDimVsTimeGPFA(seqTrain, DataSetList(nData).params);
%     setPrint(4*6, ceil(xDim/4)*5, [PlotDir 'GPFA_' DataSetList(nData).name], 'tif')
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12.3 GPFA for collected data 3d traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% runIdx               = 14;
% 
% for nData            = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     dat              = getGPFAData(nDataSet, totTargets);
%     runIdx           = runIdx + 1;
%     xDim             = ceil(size(dat(1).spikes,1)/3*2);
%     result           = neuralTraj(runIdx, dat, 'method', 'gpfa', 'xDim', xDim);
%     [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
%     plot3DGPFA(seqTrain, DataSetList(nData).params);
%     savefig([PlotDir 'GPFA3DTrace_' DataSetList(nData).name]);
% end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 12.4 GPFA for collected data EV and time constant
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% runIdx               = 14;
% 
% for nData            = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     dat              = getGPFAData(nDataSet, totTargets);
%     runIdx           = runIdx + 1;
%     xDim             = ceil(size(dat(1).spikes,1)/3*2);
%     result           = neuralTraj(runIdx, dat, 'method', 'gpfa', 'xDim', xDim);
%     [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
%     method           = plotPredErrorVsDim(runIdx, kernSD,'plotOn',false);
%     timeConstant     = DataSetList(nData).params.binsize ./ sqrt(estParams.gamma);
%     figure;
%     [ax,hBar,hLine] = plotyy(1:5,diff([0 method(5).expVar(1:5)])*100,1:5,timeConstant(1:5),'bar','plot');
%     xlabel('Latent Dimension','fontsize',12)
%     ylabel(ax(1),'% Expained Variance')%,'fontsize',12)
%     ylabel(ax(2),'Time Constant (s)')%,'fontsize',12);
%     set(ax(1),'Ylim',[-5 15]);
%     set(ax(2),'Ylim',[0  2.5]);
%     set(ax(1),'Ytick',-5:5:15);
%     set(ax(2),'Ytick',0:0.5:2.5);
%     set(ax(2),'YColor','r');
%     hBar.FaceColor  = 'k';
%     hLine.LineWidth = 2;
%     hLine.Color = [1,0,0];
%     box off
%     setPrint(6, 4, [PlotDir 'GPFAEV_' DataSetList(nData).name], 'pdf')
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12.4 GPFA for collected data EV and time constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
runIdx               = 19;

for nData            = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    dat              = getGPFAData(nDataSet, totTargets);
    runIdx           = runIdx + 1;
    xDim             = 6;
    neuralTraj(runIdx, dat, 'method', 'gpfa', 'xDim', xDim, 'numFolds', numFolds);
    result           = neuralTraj(runIdx, dat, 'method', 'gpfa', 'xDim', xDim);
    [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
    method           = plotPredErrorVsDim(runIdx, kernSD,'plotOn',false);
    timeConstant     = DataSetList(nData).params.binsize ./ sqrt(estParams.gamma);
    figure;
    [ax,hBar,hLine] = plotyy(1:xDim,diff([0 method(5).expVar(1:xDim)])*100,1:xDim,timeConstant(1:xDim),'bar','plot');
    xlabel('Latent Dimension','fontsize',12)
    ylabel(ax(1),'% Expained Variance')%,'fontsize',12)
    ylabel(ax(2),'Time Constant (s)')%,'fontsize',12);
    set(ax(1),'Ylim',[-5 15]);
    set(ax(2),'Ylim',[0  2.5]);
    set(ax(1),'Ytick',-5:5:15);
    set(ax(2),'Ytick',0:0.5:2.5);
    set(ax(2),'YColor','r');
    hBar.FaceColor  = 'k';
    hLine.LineWidth = 2;
    hLine.Color = [1,0,0];
    box off
    setPrint(6, 4, [PlotDir 'GPFAEV_' DataSetList(nData).name], 'pdf')
end