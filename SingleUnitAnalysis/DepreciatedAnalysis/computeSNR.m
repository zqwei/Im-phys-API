%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);


for nData             = [1 3 4]
    load([TempDatDir DataSetList(nData).name '.mat'])
    matSNR = nan(length(nDataSet), 2);
    for nUnit = 1:length(nDataSet)
        matSNR(nUnit, 1) = var(mean(nDataSet(nUnit).unit_yes_trial))/mean(var(nDataSet(nUnit).unit_yes_trial));
        matSNR(nUnit, 2) = var(mean(nDataSet(nUnit).unit_no_trial))/mean(var(nDataSet(nUnit).unit_no_trial));
    end
    figure;
    hold on
    [f, xi] = hist(matSNR(:, 1),100);
    semilogx(xi, f, '-b')
    [f, xi] = hist(matSNR(:, 2),100);
    semilogx(xi, f, '-r')
    hold off
    xlim([0 0.5])
   
end
