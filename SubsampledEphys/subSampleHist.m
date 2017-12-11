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

empY        = hist(mean(frData, 2), 0:1:200); % mean at 5.254 Hz
empY        = empY/numNeuron;

% low firing rate 1Hz

for frThres = [1 4 10] * 2.6 % spike count in this case
    Y           = poisspdf(0:1:200,frThres);
    Y           = Y/sum(Y);
    empRand     = rand(numNeuron, 1);
    subSampleratio = Y./empY;
    empFR       = ceil(mean(frData, 2));
    validData   = empFR < 40+frThres;
    empRatio    = subSampleratio(empFR);
    countsNum   = hist(mean(frData(empRand < empRatio' & validData, :), 2), 0:1:200);
    hold on
    stairs((0:1:200)/2.6, countsNum/sum(countsNum), 'linewid', 2.0)
    disp(mean(mean(frData(empRand < empRatio' & validData, :), 2))/2.6);  
end

xlim([0 25])

