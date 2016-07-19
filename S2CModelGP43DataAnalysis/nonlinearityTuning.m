addpath('../Func');
setDir;

numTrials        = 50;

load([TempDatDir 'Modeled_Ca_GP43_No_Noise.mat']);
s2cDataSet = nDataSet;
s2cFiringRates   = generateDPCADataLinearData(s2cDataSet, numTrials);
nls2cFiringRates = s2cFiringRates;
s2cFiringRates   = reshape(s2cFiringRates, size(s2cFiringRates, 1), []);


load([TempDatDir 'Shuffle_Ca_Slow_Short_Delay_withOLRemoval.mat']);
caDataSet        = nDataSet;
caFiringRates    = generateDPCAData(caDataSet, numTrials);
targetRates      = caFiringRates;
caFiringRates    = reshape(caFiringRates, size(caFiringRates, 1), []);

rhoMat                 = corr(caFiringRates', s2cFiringRates', 'type', 'Spearman');
[~, maxIndex]   = max(rhoMat, [], 1);

f = @(p,x) p(1) + p(2)./ (1 + 10.^((p(3)-x)*p(4)));
g = @(p,x) p(1) + p(2)./ (1 + (p(3)./x).^p(4));

nlParams               = nan(size(s2cFiringRates, 1), 4);

for nUnit = 1:size(s2cFiringRates, 1)
    xdata = s2cFiringRates(nUnit, :);
    ydata = caFiringRates(maxIndex(nUnit), :);
    [param,~]=sigm_fit(xdata, ydata, [], [], false);
    param(2) = param(2) - param(1);
%     % remove nlinfit, which generate complex number in fit
% %     param   = nlinfit(xdata,ydata,g,[min(ydata) max(ydata)-min(ydata) (max(ydata)+min(ydata))/2 1]);
%     param   = lsqcurvefit(g, [min(ydata) max(ydata)-min(ydata) (max(ydata)+min(ydata))/2 1], xdata, ydata);
    nlParams(nUnit, :) = param;
%     nls2cFiringRates(nUnit, :, :, :) = g(param, nls2cFiringRates(nUnit, :, :, :));
    nls2cFiringRates(nUnit, :, :, :) = f(param, nls2cFiringRates(nUnit, :, :, :));
end

noiseRates{1} = targetRates(maxIndex, :, :, :) - nls2cFiringRates;



% load([TempDatDir 'Modeled_Ca_Short_Decay_No_Noise.mat']);
% s2cDataSet = nDataSet;
% s2cFiringRates   = generateDPCADataLinearData(s2cDataSet, numTrials);
% nls2cFiringRates = f(param, s2cFiringRates);
% s2cFiringRates   = reshape(nls2cFiringRates, size(s2cFiringRates, 1), []);
% 
% load([TempDatDir 'Shuffle_Ca_Slow_Short_Delay_withOLRemoval.mat']);
% caDataSet        = nDataSet;
% caFiringRates    = generateDPCAData(caDataSet, numTrials);
% targetRates      = caFiringRates;
% caFiringRates    = reshape(caFiringRates, size(caFiringRates, 1), []);
% rhoMat                 = corr(caFiringRates', s2cFiringRates', 'type', 'Spearman');
% [~, maxIndex]   = max(rhoMat, [], 1);
% noiseRates{2} = targetRates(maxIndex, :, :, :) - nls2cFiringRates;

save('FineTunedNLParams', 'noiseRates', 'nlParams');



