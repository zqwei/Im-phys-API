%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.  Collected population decision decodability over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.  Collected population decision decodability over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_6_1
% The same amount of units in analysis
addpath('../Func');
setDir;

ROCThres            = 0.6;

for nData             = [3 6];
    load([TempDatDir DataSetList(nData).name '.mat'])
    figure;
    selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';% & [DataSetList(nData).cellinfo(:).cellType] == 1;
    selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = nDataSet(selectedNeuronalIndex);
    
    m                     = min(ceil(sqrt(length(nDataSet))),10);
    
    for nNeuron           = 1:min(length(nDataSet),100)
        mean_yes      = mean(nDataSet(nNeuron).unit_yes_trial);
        mean_no       = mean(nDataSet(nNeuron).unit_no_trial);
        var_yes       = sem(nDataSet(nNeuron).unit_yes_trial);
        var_no        = sem(nDataSet(nNeuron).unit_no_trial);
        subplot(m, m, nNeuron)
        hold on;
        shadedErrorBar(params.timeSeries, mean_yes, var_yes, {'-b', 'linewid', 1.0}, 0.5);
        shadedErrorBar(params.timeSeries, mean_no, var_no, {'-r', 'linewid', 1.0}, 0.5);
        xlim([params.timeSeries(1) params.polein]);
%         % ylim([yAxes_set(1) yAxes_set(2)]);
%         gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
        hold off; 
        xlabel('Time (s)');
        ylabel('DF/F')
    end
    

end