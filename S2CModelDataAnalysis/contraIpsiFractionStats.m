%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% no-selective cell      : no selectivity -- 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffleConfounding.mat']);
DataSetListAnm = DataSetList([1 2 3 6]); 
load ([TempDatDir 'DataListShuffle.mat']);

for nData = [1 3 4]
    if nData == 1
        load([TempDatDir DataSetList(nData).name '.mat'])
        depth = [DataSetList(nData).cellinfo.depth];
        validDepth = depth<700 & depth>100;
        nDataSet   = nDataSet(validDepth);
        [~, ~, anmIndex] = unique(cell2mat({DataSetListAnm(nData).cellinfo(validDepth).anmName}'), 'rows');
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
        [~, ~, anmIndex] = unique(cell2mat({DataSetListAnm(nData).cellinfo(~neuronRemoveList).anmName}'), 'rows');
    end
    logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
    unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);
    params      = DataSetList(nData).params;
    contraIndex = false(length(nDataSet), 1);
    cellType    = [DataSetList(nData).cellinfo.cellType]';
    yesActMat   = nan(length(nDataSet), length(params.timeSeries));
    noActMat    = nan(length(nDataSet), length(params.timeSeries));
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);

    for nUnit   = 1:length(nDataSet)
        yesTrial = mean(nDataSet(nUnit).unit_yes_trial);
        noTrial  = mean(nDataSet(nUnit).unit_no_trial);
        yesActMat(nUnit, :)  = yesTrial;
        noActMat(nUnit, :)   = noTrial;
        contraIndex(nUnit)   = sum(noTrial(timePoints(2):end))<sum(yesTrial(timePoints(2):end));
    end
    
    contraIndexANM = nan(anmIndex(end),1);
    
    for nAnm    = 1:anmIndex(end)
        if sum(anmIndex == nAnm) > 50
            contraIndexANM(nAnm) = mean(contraIndex & unitGroup>0 & anmIndex == nAnm)/mean(unitGroup>0 & anmIndex == nAnm);
        end
    end
    
    contraIndexData {nData} = contraIndexANM;
    meanContraIndexData(nData) = mean(contraIndex & unitGroup>0)/mean(unitGroup>0);
end

nData = 1;
load([TempDatDir DataSetList(nData).name '.mat'])
depth = [DataSetList(nData).cellinfo.depth];
validDepth = depth<700 & depth>100;
nDataSet   = nDataSet(validDepth);
[~, ~, anmIndex] = unique(cell2mat({DataSetListAnm(nData).cellinfo(validDepth).anmName}'), 'rows');
load ([TempDatDir 'DataListS2CModel.mat']);
for nData = 1:2
    load([TempDatDir DataSetList(nData).name '.mat'])
    validDepth = depth<700 & depth>100;
    nDataSet   = nDataSet(validDepth);
    logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
    unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);
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
        contraIndex(nUnit)   = sum(noTrial(timePoints(2):end))<sum(yesTrial(timePoints(2):end));
    end
    
    contraIndexANM = nan(anmIndex(end),1);
    
    for nAnm    = 1:anmIndex(end)
        if sum(anmIndex == nAnm) > 50
            contraIndexANM(nAnm) = mean(contraIndex & unitGroup>0 & anmIndex == nAnm)/mean(unitGroup>0 & anmIndex == nAnm);
        end
    end
    
    contraIndexData {nData+4} = contraIndexANM;
    meanContraIndexData(nData+4) = mean(contraIndex & unitGroup>0)/mean(unitGroup>0);
end

meanContraIndexData(2) = [];
contraIndexData(2) = [];
stdContraIndexData = cellfun(@nanstd, contraIndexData)./sqrt(cellfun(@(x) sum(~isnan(x)), contraIndexData));


figure('visible', 'on');
hold on
barerror((1:5)', meanContraIndexData', stdContraIndexData', 1, ['k'],  ['k'])
gridxy([], [0.5], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
xlim([0.5 5.5])
ylim([0.4 0.65])
xlabel('Dataset index')
ylabel('Frac. Contra')
box off;
setPrint(8, 6, [PlotDir 'S2CModel/ContraFraction'])