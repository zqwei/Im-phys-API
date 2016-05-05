%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% no-selective cell      : no selectivity -- 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CModel.mat']);
load('spikeContraIndex.mat')

for nData = [3 4]
    load([TempDatDir DataSetList(nData).name '.mat'])
    logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
    unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);
    params      = DataSetList(nData).params;
%     contraIndex = false(length(nDataSet), 1);
    yesActMat   = nan(length(nDataSet), length(params.timeSeries));
    noActMat    = nan(length(nDataSet), length(params.timeSeries));
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    contraIndex = spikeContraIndex;

    for nUnit   = 1:length(nDataSet)
        yesTrial = mean(nDataSet(nUnit).unit_yes_trial);
        noTrial  = mean(nDataSet(nUnit).unit_no_trial);
        yesActMat(nUnit, :)  = yesTrial;
        noActMat(nUnit, :)   = noTrial;
%         contraIndex(nUnit)   = sum(noTrial(timePoints(2):end))>sum(yesTrial(timePoints(2):end));
    end    
%     mean(contraIndex)
    
    figure;
    subplot(1, 2, 1)
    hold on
    shadedErrorBar(params.timeSeries, mean(noActMat(contraIndex,:)), sem(noActMat(contraIndex,:)),'-b')
    shadedErrorBar(params.timeSeries, mean(yesActMat(contraIndex,:)), sem(yesActMat(contraIndex,:)),'-r')
    title('contra neuron only')
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    xlim([params.timeSeries(1) params.timeSeries(end)])
    xlabel('Time (ms)')
    box off
    ylabel('Mean activity')

    subplot(1, 2, 2)
    hold on
    shadedErrorBar(params.timeSeries, mean(noActMat(~contraIndex,:)), sem(noActMat(~contraIndex,:)),'-b')
    shadedErrorBar(params.timeSeries, mean(yesActMat(~contraIndex,:)), sem(yesActMat(~contraIndex,:)),'-r')
    title('ipsi neuron only')
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    xlim([params.timeSeries(1) params.timeSeries(end)])
    xlabel('Time (ms)')
    box off
    ylabel('Mean activity')
    
    setPrint(8*2, 6, [PlotDir 'S2CModel\SingleUnitsContraIpsi_' DataSetList(nData).name])
    
    
    for nUnit   = 1:length(nDataSet)
        yesTrial = mean(nDataSet(nUnit).unit_yes_trial_linear);
        noTrial  = mean(nDataSet(nUnit).unit_no_trial_linear);
        yesActMat(nUnit, :)  = yesTrial;
        noActMat(nUnit, :)   = noTrial;
%         contraIndex(nUnit)   = sum(noTrial(timePoints(2):end))>sum(yesTrial(timePoints(2):end));
    end    
%     mean(contraIndex)
    
    figure;
    subplot(1, 2, 1)
    hold on
    shadedErrorBar(params.timeSeries, mean(noActMat(contraIndex,:)), sem(noActMat(contraIndex,:)),'-b')
    shadedErrorBar(params.timeSeries, mean(yesActMat(contraIndex,:)), sem(yesActMat(contraIndex,:)),'-r')
    title('contra neuron only')
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    xlim([params.timeSeries(1) params.timeSeries(end)])
    xlabel('Time (ms)')
    box off
    ylabel('Mean activity')

    subplot(1, 2, 2)
    hold on
    shadedErrorBar(params.timeSeries, mean(noActMat(~contraIndex,:)), sem(noActMat(~contraIndex,:)),'-b')
    shadedErrorBar(params.timeSeries, mean(yesActMat(~contraIndex,:)), sem(yesActMat(~contraIndex,:)),'-r')
    title('ipsi neuron only')
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    xlim([params.timeSeries(1) params.timeSeries(end)])
    xlabel('Time (ms)')
    box off
    ylabel('Mean activity')
    
    setPrint(8*2, 6, [PlotDir 'S2CModel\SingleUnitsContraIpsiLinear_' DataSetList(nData).name])
    
end
% figure;
% hold on
% shadedErrorBar(params.timeSeries, mean(noActMat(contraIndex & cellType == 1 & unitGroup>0,:)), sem(noActMat(contraIndex & cellType == 1 & unitGroup>0,:)),'-b')
% shadedErrorBar(params.timeSeries, mean(yesActMat(contraIndex & cellType == 1 & unitGroup>0,:)), sem(yesActMat(contraIndex & cellType == 1 & unitGroup>0,:)),'-r')
% title('contra selective pyr neuron')
% 
% figure;
% hold on
% shadedErrorBar(params.timeSeries, mean(noActMat(~contraIndex & cellType == 1 & unitGroup>0,:)), sem(noActMat(~contraIndex & cellType == 1 & unitGroup>0,:)),'-b')
% shadedErrorBar(params.timeSeries, mean(yesActMat(~contraIndex & cellType == 1 & unitGroup>0,:)), sem(yesActMat(~contraIndex & cellType == 1 & unitGroup>0,:)),'-r')
% title('ipsi selective pyr neuron')
% 
% figure;
% hold on
% shadedErrorBar(params.timeSeries, mean(noActMat(contraIndex & cellType == 0 & unitGroup>0,:)), sem(noActMat(contraIndex & cellType == 0 & unitGroup>0,:)),'-b')
% shadedErrorBar(params.timeSeries, mean(yesActMat(contraIndex & cellType == 0 & unitGroup>0,:)), sem(yesActMat(contraIndex & cellType == 0 & unitGroup>0,:)),'-r')
% title('contra selective FS neuron')
% 
% figure;
% hold on
% shadedErrorBar(params.timeSeries, mean(noActMat(~contraIndex & cellType == 0 & unitGroup>0,:)), sem(noActMat(~contraIndex & cellType == 0 & unitGroup>0,:)),'-b')
% shadedErrorBar(params.timeSeries, mean(yesActMat(~contraIndex & cellType == 0 & unitGroup>0,:)), sem(yesActMat(~contraIndex & cellType == 0 & unitGroup>0,:)),'-r')
% title('ipsi selective FS neuron')
% 
% figure;
% hold on
% shadedErrorBar(params.timeSeries, mean(noActMat(contraIndex & unitGroup>0,:)), sem(noActMat(contraIndex & unitGroup>0,:)),'-b')
% shadedErrorBar(params.timeSeries, mean(yesActMat(contraIndex & unitGroup>0,:)), sem(yesActMat(contraIndex & unitGroup>0,:)),'-r')
% title('contra selective neuron')
% 
% figure;
% hold on
% shadedErrorBar(params.timeSeries, mean(noActMat(~contraIndex & unitGroup>0,:)), sem(noActMat(~contraIndex & unitGroup>0,:)),'-b')
% shadedErrorBar(params.timeSeries, mean(yesActMat(~contraIndex & unitGroup>0,:)), sem(yesActMat(~contraIndex & unitGroup>0,:)),'-r')
% title('ipsi selective neuron')



% figure;
% hold on
% shadedErrorBar(params.timeSeries, mean(noActMat), sem(noActMat),'-b')
% shadedErrorBar(params.timeSeries, mean(yesActMat), sem(yesActMat),'-r')
% title('all neuron')
% 
% figure;
% hold on
% shadedErrorBar(params.timeSeries, mean(noActMat(unitGroup>0,:)), sem(noActMat(unitGroup>0,:)),'-b')
% shadedErrorBar(params.timeSeries, mean(yesActMat(unitGroup>0,:)), sem(yesActMat(unitGroup>0,:)),'-r')
% title('selective neuron only')
% 
% figure;
% hold on
% shadedErrorBar(params.timeSeries, mean(noActMat(cellType == 0,:)), sem(noActMat(cellType == 0,:)),'-b')
% shadedErrorBar(params.timeSeries, mean(yesActMat(cellType == 0,:)), sem(yesActMat(cellType == 0,:)),'-r')
% title('FS neuron only')
% 
% figure;
% hold on
% shadedErrorBar(params.timeSeries, mean(noActMat(cellType == 1,:)), sem(noActMat(cellType == 1,:)),'-b')
% shadedErrorBar(params.timeSeries, mean(yesActMat(cellType == 1,:)), sem(yesActMat(cellType == 1,:)),'-r')
% title('Pyr neuron only')
