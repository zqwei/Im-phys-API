%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating new data using fine tuned params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;

% measure of nonlinearity
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');

dataSetNames  = {'6sAAV', 'GP43'};

for nData     = 1:2
    nlParamsnData = squeeze(nlParams(nData, :, :));
    mdlr          = fitlm(log10(nlParamsnData(:,3)), log10(nlParamsnData(:, 4)),'RobustOpts','on');
    p1            = mdlr.Coefficients.Estimate(2);
    p0            = mdlr.Coefficients.Estimate(1);
    x_est         = linspace(min(log10(nlParamsnData(:,3))), max(log10(nlParamsnData(:,3))), 101);
    y_est         = 10.^(x_est * p1 + p0);
    x_est         = 10.^(x_est);
    figure;
    loglog(nlParams_1(:, 3), nlParams_1(:, 4), 'o')
    hold on
    loglog(x_est, y_est)
    box off
    xlabel('EC50')
    ylabel('NL')
    xlim([min(nlParamsnData(:,3)), max(nlParamsnData(:,3))])
    ylim([min(nlParamsnData(:,4)), max(nlParamsnData(:,4))])
    setPrint(8, 6, ['NLParams_' dataSetNames{nData}], 'pdf')
    disp([p0, p1, max(nlParamsnData(:,3))]);
end