%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% no-selective cell      : no selectivity -- 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load([TempDatDir 'DataListShuffle.mat']);

figure;
hold on
nIndex = 0;

for nData = [1 10]
    if nData == 1
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end
    depth_list          = [nDataSet.depth_in_um]';    
    oldDataSet          = nDataSet;
    selectedNeuronalIndex = depth_list < 471;
    nDataSet              = oldDataSet(selectedNeuronalIndex);
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params); 
    params      = DataSetList(nData).params;
    contraIndex = false(length(nDataSet), 1);
    yesActMat   = nan(length(nDataSet), length(params.timeSeries));
    noActMat    = nan(length(nDataSet), length(params.timeSeries));
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);

    for nUnit   = 1:length(nDataSet)
        yesTrial = mean(nDataSet(nUnit).unit_yes_trial);
        noTrial  = mean(nDataSet(nUnit).unit_no_trial);
        yesActMat(nUnit, :)  = yesTrial;
        noActMat(nUnit, :)   = noTrial;
        if nData ==1 
            contraIndex(nUnit)   = sum(noTrial(timePoints(2):timePoints(5)))<sum(yesTrial(timePoints(2):timePoints(5)));
        else
            contraIndex(nUnit)   = sum(noTrial(timePoints(2):timePoints(4)))<sum(yesTrial(timePoints(2):timePoints(4)));    
        end
    end
    
    nIndex = nIndex + 1;
    [~, ~, anmIndex] = unique(cell2mat({DataSetList(nData).cellinfo.anmName}'), 'rows');
    anmIndex  = anmIndex(~neuronRemoveList);
    anmIndex = anmIndex(selectedNeuronalIndex);
    if nData == 1
        anmSpkIndex = anmIndex;
    end
    
    contraCount = grpstats(unitGroup~=0 & contraIndex, anmIndex, 'sum');
    ipsiCount = grpstats(unitGroup~=0 & ~contraIndex, anmIndex, 'sum');
    
    bar(nIndex, sum(contraCount)/(sum(contraCount)+sum(ipsiCount)), ...
        'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    errorbar(nIndex, sum(contraCount)/(sum(contraCount)+sum(ipsiCount)), ...
        sem(contraCount./(contraCount+ipsiCount)), 'k')
end

load([TempDatDir 'DataListS2C6fModel.mat']);

for nData = 2
    load([TempDatDir DataSetList(nData).name '.mat'])
        
%     logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
%     unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params); 
    params      = DataSetList(nData).params;
    contraIndex = false(length(nDataSet), 1);
    yesActMat   = nan(length(nDataSet), length(params.timeSeries));
    noActMat    = nan(length(nDataSet), length(params.timeSeries));
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);

    for nUnit   = 1:length(nDataSet)
        yesTrial = mean(nDataSet(nUnit).unit_yes_trial);
        noTrial  = mean(nDataSet(nUnit).unit_no_trial);
        yesActMat(nUnit, :)  = yesTrial;
        noActMat(nUnit, :)   = noTrial;
        contraIndex(nUnit)   = sum(noTrial(timePoints(2):timePoints(5)))<sum(yesTrial(timePoints(2):timePoints(5)));
    end
    
    nIndex = nIndex + 1;
    contraCount = grpstats(unitGroup~=0 & contraIndex, anmSpkIndex, 'sum');
    ipsiCount = grpstats(unitGroup~=0 & ~contraIndex, anmSpkIndex, 'sum');
    
    bar(nIndex, sum(contraCount)/(sum(contraCount)+sum(ipsiCount)), ...
        'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    errorbar(nIndex, sum(contraCount)/(sum(contraCount)+sum(ipsiCount)), ...
        sem(contraCount./(contraCount+ipsiCount)), 'k')
    
end

set(gca, 'TickDir', 'out')
ylim([0.4 0.64])
xlim([0.5 nIndex+0.5])
setPrint(8, 10, [PlotDir 'SingleUnitsContraIpsi\Faction'])