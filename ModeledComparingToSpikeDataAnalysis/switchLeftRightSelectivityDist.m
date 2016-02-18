%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% left/right cell      : change of selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);


if ~exist([PlotDir 'ModeledComparingSingleUnitsTscore'],'dir')
    mkdir([PlotDir 'ModeledComparingSingleUnitsTscore'])
end

nData          = 1;
load([TempDatDir DataSetList(nData).name '.mat'])

spkPre         = nan(length(nDataSet), 1);
spkLeft        = nan(length(nDataSet), 1);
spkRight       = nan(length(nDataSet), 1);

params         = DataSetList(nData).params;
timePoints     = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);


for nNeuron    = 1:length(nDataSet)
    nNeuronData        = nDataSet(nNeuron);
    spkPre(nNeuron)    = mean(mean([nNeuronData.unit_yes_trial(:, timePoints(1):timePoints(2)); nNeuronData.unit_no_trial(:, timePoints(1):timePoints(2))]));
    spkLeft(nNeuron)   = mean(mean(nNeuronData.unit_yes_trial(:, timePoints(2):timePoints(end))));
    spkRight(nNeuron)  = mean(mean(nNeuronData.unit_no_trial(:, timePoints(2):timePoints(end))));
end

isLeft  = spkLeft > spkRight;
isRight = spkLeft < spkRight;
actMat  = abs(max(spkLeft, spkRight) - spkPre);


nGroup = round(actMat);

[probLeft, ~, groupIndexLeft]= grpstats(isLeft, nGroup, {'sum', 'numel', 'gname'});
groupIndexLeft = str2double(groupIndexLeft);
[probRight, ~, groupIndexRight]= grpstats(isRight, nGroup, {'sum', 'numel', 'gname'});
groupIndexRight = str2double(groupIndexRight);

figure;
subplot(1, 2, 1)
plot(groupIndexLeft, probLeft, '-o', 'linewid', 2)

subplot(1,2, 2)
plot(groupIndexRight, probRight, '-o', 'linewid', 2)


for nData = 2:length(DataSetList)
    
    load([TempDatDir DataSetList(nData).name '.mat']);
    for nNeuron    = 1:length(nDataSet)
        nNeuronData        = nDataSet(nNeuron);
        spkLeft(nNeuron)   = mean(mean(nNeuronData.unit_yes_trial(:, timePoints(2):timePoints(end))));
        spkRight(nNeuron)  = mean(mean(nNeuronData.unit_no_trial(:, timePoints(2):timePoints(end))));
    end
    
    isLeft = spkLeft > spkRight;
    isRight = spkLeft < spkRight;
    [probLeft, ~, groupIndexLeft]= grpstats(isLeft, nGroup, {'sum', 'numel', 'gname'});
    groupIndexLeft = str2double(groupIndexLeft);
    [probRight, ~, groupIndexRight]= grpstats(isRight, nGroup, {'sum', 'numel', 'gname'});
    groupIndexRight = str2double(groupIndexRight);

    
    subplot(1, 2, 1);hold on
    plot(groupIndexLeft, probLeft, '-o', 'linewid', 2)

    subplot(1,2, 2);hold on
    plot(groupIndexRight, probRight, '-o', 'linewid', 2)
end

subplot(1, 2, 1)
xlabel('|R_{max} - R_{pre}| (Hz)')
ylabel('# left-seletive neurons')
% xlim([-0.1 15])
% ylim([0 1])
set(gca, 'TickDir', 'out')
box off

subplot(1, 2, 2)
xlabel('|R_{max} - R_{pre}| (Hz)')
ylabel('# right-seletive neurons')
% xlim([-0.1 15])
% ylim([0 1])
set(gca, 'TickDir', 'out')
box off

legend({'Spike','Ca long decay', 'Ca short decay'})

setPrint(8*2, 6, 'LeftRightNeuronNumber', 'tif')