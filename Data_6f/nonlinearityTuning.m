addpath('../Func');
setDir;

numTrials        = 50;
load([TempDatDir 'Modeled_GP517.mat']);
s2cDataSet = nDataSet;
s2cFiringRates   = generateDPCADataLinearData(s2cDataSet, numTrials);
s2cFiringRatesMean   = squeeze(mean(s2cFiringRates, 4));

load([TempDatDir 'Shuffle_Ca_Fast_SShort_Delay_withOLRemoval.mat']);
caDataSet        = nDataSet;
caFiringRates    = generateDPCAData(caDataSet, numTrials);
caFiringRatesIpsi= caFiringRates;
caFiringRatesIpsi(:, 1, :, :) = caFiringRates(:, 2, :, :);
caFiringRatesIpsi(:, 2, :, :) = caFiringRates(:, 1, :, :);
caFiringRatesMean = squeeze(mean(caFiringRates, 4));
caFiringRatesIpsiMean = squeeze(mean(caFiringRatesIpsi, 4));
caFiringRatesMean    = [squeeze(caFiringRatesMean(:, 1, :)), squeeze(caFiringRatesMean(:, 2, :))];
numCaTime            = size(caFiringRatesMean, 2)/2;
caFiringRatesIpsiMean= [squeeze(caFiringRatesIpsiMean(:, 1, :)), squeeze(caFiringRatesIpsiMean(:, 2, :))];
s2cFiringRatesMean   = [squeeze(s2cFiringRatesMean(:, 1, 1:numCaTime)), squeeze(s2cFiringRatesMean(:, 2, 1:numCaTime))];

rhoMat                 = corr(caFiringRatesMean', s2cFiringRatesMean', 'type', 'Spearman');
rhoMatIpsi             = corr(caFiringRatesIpsiMean', s2cFiringRatesMean', 'type', 'Spearman');
[maxRho, maxIndex]   = max(rhoMat, [], 1);
[maxRhoIpsi, maxIndexIpsi]   = max(rhoMatIpsi, [], 1);

f = @(p,x) p(1) + p(2)./ (1 + exp((p(3)-x)*p(4)));
nlParams = nan(size(s2cFiringRatesMean, 1), 4);
for nUnit = 1:size(s2cFiringRatesMean, 1)
    xdata = s2cFiringRates(nUnit, :, 1:numCaTime, :);
    xdata = sort(xdata, 4);
    if maxRho(nUnit)>maxRhoIpsi(nUnit)
        ydata = caFiringRates(maxIndex(nUnit), :, :, :);
    else
        ydata = caFiringRatesIpsi(maxIndexIpsi(nUnit), :, :, :);
    end
    ydata = sort(ydata, 4);
    init_para = [min(ydata(:)), max(ydata(:)), median(xdata(:)), xdata(:)\ydata(:)];
    param = lsqcurvefit(f,init_para,xdata(:),ydata(:),...
            [min(ydata(:)), max(ydata(:))/2, median(xdata(:))/3, 0], ...
            [0, max(ydata(:)), median(xdata(:))*3, 20]);
    if param(4)<=0 || param(3)<0
        param = nan(1, 4);
    end
    nlParams(nUnit, :) = param;
end


f = @(p,x) p(1) + p(2)./ (1 + exp((p(3)-x)*p(4)));
nlParams = nan(size(s2cFiringRatesMean, 1), 4);
for nUnit = 1:size(s2cFiringRatesMean, 1)
    xdata = s2cFiringRatesMean(nUnit, :);
    if maxRho(nUnit)>maxRhoIpsi(nUnit)
        ydata = caFiringRatesMean(maxIndex(nUnit), :);
    else
        ydata = caFiringRatesIpsiMean(maxIndexIpsi(nUnit), :);
    end
    [param,~]=sigm_fit(xdata, ydata, [], [], false);
    param(2) = param(2) - param(1);
    param(4) = param(4) * log(10);
    if param(4)<=0 || param(3)<0
        param = nan(1, 4);
    end
    nlParams(nUnit, :) = param;
    nls2cFiringRates(nUnit, :, :, :) = f(param, s2cFiringRates(nUnit, :, :, :));
end


% examine fitted traces
for nUnit = 1:size(s2cFiringRates, 1)
    xdata = nls2cFiringRates(nUnit, :, :, :);
    xdata = squeeze(mean(xdata, 4));
    figure('visible', 'on');
    subplot(1, 2, 1)
    plot(xdata')
    box off
    rhoNCell = max(maxRho(nUnit), maxRhoIpsi(nUnit));
    if maxRho(nUnit)>maxRhoIpsi(nUnit)
        ydata = caFiringRatesMean(maxIndex(nUnit), :);
    else
        ydata = caFiringRatesIpsiMean(maxIndexIpsi(nUnit), :);
    end
    subplot(1, 2, 2)
    plot([ydata(1:numCaTime)', ydata(numCaTime+1:numCaTime*2)'])
    box off
    title(num2str(rhoNCell))
%     setPrint(8*2, 6, ['NLFitPlots/Cell_id_' num2str(nUnit, '%03d')])
%     close all
end

save([TempDatDir 'FineTuned6fNLParams.mat'], 'nlParams');