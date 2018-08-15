% Examine the firing rates in ephys vs whole cell recording data

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
load([TempDatDir DataSetList(1).name '.mat'])

% ephys
timePeriod  = 8:47;
numNeuron   = length(nDataSet);

figure;
frData      = nan(numNeuron, 2);

for nUnit   = 1:numNeuron
    frData(nUnit, 1) = mean(mean(nDataSet(nUnit).unit_yes_trial(:, timePeriod)))*2.6;
    frData(nUnit, 2) = mean(mean(nDataSet(nUnit).unit_no_trial(:, timePeriod)))*2.6;
end

frData      = mean(frData, 2);

empY        = hist(frData, 0:1:200); % mean at 5.254 Hz
empY        = empY/numNeuron;

stairs((0:1:200)/2.6, empY, 'linewid', 1.0, 'color', 'k')

numFold     = 30;

% low firing rate 1Hz

for frThres = [1 4 10] % spike count in this case
    Y           = poisspdf(0:1:200,frThres * 2.6 );
    Y           = Y/sum(Y);
    subSampleratio = Y./empY;
    empFR       = ceil(frData);
    validData   = empFR < 40+frThres * 2.6 ;
    empRatio    = subSampleratio(empFR);
    totCounts   = 0;
    meanFR      = 0;
    
    validMat    = false(numNeuron, numFold);
%     load(['validMat_' num2str(frThres, '%02d')], 'validMat')
    for nFold   = 1:numFold
        empRand     = rand(numNeuron, 1);
        validMat(:, nFold) = empRand < empRatio' & validData;
        countsNum   = hist(frData(validMat(:, nFold)), 0:1:200);
        countsNum   = countsNum/sum(countsNum);
        totCounts   = totCounts + countsNum/numFold;
        meanFR      = meanFR + mean(frData(validMat(:, nFold)))/numFold;
    end
    
    hold on
    stairs((0:1:200)/2.6, totCounts, 'linewid', 1.0)
    disp(meanFR/2.6); 
%     save(['validMat_' num2str(frThres, '%02d')], 'validMat')
end

xlim([0 25])
xlabel('Firing rate (/s)')
ylabel('Frac. neuron')
set(gca, 'TickDir', 'out')
setPrint(8, 6, 'Firing_rate_distribution', 'pdf')
