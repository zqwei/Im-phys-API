%
% Compute variance - signal correlation in ephys data
% 
% -------------------------------------------------------------------------
% version 1.0
%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
load([TempDatDir DataSetList(1).name '.mat'])
params    = DataSetList(1).params;
binSize   = params.binsize;

numTime   = length(params.timeSeries);
numNeuron = length(nDataSet);
timePoints= timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);

% instanteous noise level
var_mat   = nan(numNeuron, numTime * 2); 
sig_mat   = nan(numNeuron, numTime * 2); 

for nNeuron = 1:numNeuron
   sig_mat(nNeuron, 1:numTime)     = mean(nDataSet(nNeuron).unit_yes_trial);
   var_mat(nNeuron, 1:numTime)     = var(nDataSet(nNeuron).unit_yes_trial);
   sig_mat(nNeuron, 1+numTime:end) = mean(nDataSet(nNeuron).unit_no_trial);
   var_mat(nNeuron, 1+numTime:end) = var(nDataSet(nNeuron).unit_no_trial);
end

% sig - var relation
% in log-log scale
loglog(sig_mat(:), sqrt(var_mat(:))./sig_mat(:), '.')
% fit
% p1 * x + p2
% p1 = 0.9482
% p2 = 2.79

% distribution of var:sig ratio at finer time points
ksdensity(var_mat(:)./sig_mat(:))
hist(var_mat(:)./sig_mat(:) * binSize,10000)
ratio_list = var_mat(:)./sig_mat(:) * binSize;
ratio_list(isnan(ratio_list)) = 1;
% ratio_list<4 nearly 99.71%
mean(ratio_list);

% same for each neuron
ratio_list = var_mat./sig_mat * binSize;
ratio_list(isnan(ratio_list)) = 1;
mean(ratio_list);

% same for each neuron -- normalized activity
% hist(mean(sig_mat* binSize, 2),50)
[f, xi] = hist(mean(sig_mat* binSize, 2),50);
func_f  = fit(xi',f','exp1');
figure; hold on
plot(xi, f, 'o')
plot(xi, feval(func_f,xi), '--')
mean(mean(sig_mat* binSize, 2))
std(mean(sig_mat* binSize, 2))

% epoch noise level
% numPeriod = length(timePoints)-2;
% var_mat   = nan(numNeuron, numPeriod * 2); 
% sig_mat   = nan(numNeuron, numPeriod * 2); 

% for nNeuron = 1:numNeuron
%     for nPeriod   = 1:numPeriod
%         actmat    = mean(nDataSet(nNeuron).unit_yes_trial(:, timePoints(nPeriod+1):timePoints(nPeriod+2)), 2);
%         sig_mat(nNeuron, nPeriod)     = mean(actmat);
%         var_mat(nNeuron, nPeriod)     = var(actmat);
%         actmat    = mean(nDataSet(nNeuron).unit_no_trial(:, timePoints(nPeriod+1):timePoints(nPeriod+2)), 2);
%         sig_mat(nNeuron, nPeriod+numPeriod) = mean(actmat);
%         var_mat(nNeuron, nPeriod+numPeriod) = var(actmat);
%     end
% end

% distribution of var:sig ratio at epoch level
% ksdensity(var_mat(:)./sig_mat(:))
% hist(var_mat(:)./sig_mat(:).^2)

% probably check this only for selective neuron
