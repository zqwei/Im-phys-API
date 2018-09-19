%
% plotMeanActivityTrace.m
% 
% based on plotMeanActivityImagesc.m
%
% ----------------------------
% Output:
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function plotMeanActivityTrace (nDataSet, nNeuron, params, ylabels, xlabels)
    mean_yes      = mean(nDataSet(nNeuron).unit_yes_trial);
    mean_no       = mean(nDataSet(nNeuron).unit_no_trial);

    var_yes       = sem(nDataSet(nNeuron).unit_yes_trial);
    var_no        = sem(nDataSet(nNeuron).unit_no_trial);
    

    hold on;
    shadedErrorBar(params.timeSeries, mean_yes, var_yes, {'-b', 'linewid', 1.0}, 0.5);
    shadedErrorBar(params.timeSeries, mean_no, var_no, {'-r', 'linewid', 1.0}, 0.5);
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off; 
    xlabel(xlabels);
    ylabel(ylabels)
    
end