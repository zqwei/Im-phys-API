%
% plotAnovaLogPSpikeExampleNeurons.m
% 
%
% Spiking dataset
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

function plotAnovaLogPSpikeExampleNeurons (logPValue, nDataSetList, nDataSet, filterInUse)

    figure;
    SpikingDataDir= '../../../Data_In_Use/Dataset_Comparison/ElectrophysiologyData/';
    
    params        = nDataSetList.params;
    thresLogP     = -log(0.05);
    numTime       = 15;
    numUnit       = size(logPValue, 1);
    minNumTrial   = 5;
    numTrialMat   = false(numUnit, 1);
    
    for nUnit     = 1:length(nDataSet)
        numTrialMat(nUnit) =    length(nDataSet(nUnit).unit_yes_trial_index) > minNumTrial && ...
                                length(nDataSet(nUnit).unit_no_trial_index ) > minNumTrial && ...
                                length(nDataSet(nUnit).unit_yes_error_index) > minNumTrial && ...
                                length(nDataSet(nUnit).unit_no_error_index ) > minNumTrial;
    end
    
    maxLogP       = 15;
    maxSelectivityMat = sum(sum(logPValue > maxLogP, 3), 2) == 0;
    numTrialMat   = numTrialMat & maxSelectivityMat; %& nDataSetList.ActiveNeuronIndex;
    timePoints    = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
%     selectivityMat = sum(logPValue(:, :, timePoints(2):end) > thresLogP + 0.5, 3) > numTime;
    unitGroup     = zeros(numUnit, 1);
    
    gAP           = sum(logPValue(:, 1, timePoints(2):timePoints(4)) > thresLogP, 3) > numTime;
    gLR           = sum(logPValue(:, 2, timePoints(2):timePoints(4)) > thresLogP, 3) > numTime;
    gCE           = sum(logPValue(:, 3, timePoints(4):timePoints(5)) > thresLogP, 3) > numTime;
    unitGroup( gAP & ~gLR & ~gCE) = 1; % group for AP
    unitGroup(~gAP &  gLR & ~gCE) = 2; % group for LR
    unitGroup(~gAP & ~gLR &  gCE) = 3; % group for CE
    unitGroup( gAP &  gLR & ~gCE) = 4; % group for AP + LR
    unitGroup( gAP & ~gLR &  gCE) = 5; % group for AP + CE
    unitGroup(~gAP &  gLR &  gCE) = 6; % group for LR + CE
    unitGroup( gAP &  gLR &  gCE) = 7; % group for AP + LR + CE 
%     unitGroup(~gAP & ~gLR & ~gCE) = 8; % group for none
    unitGroup( sum(sum(logPValue > thresLogP, 3), 2) == 0) = 8;
    
%     numUnitGroup                  = histc(unitGroup, 1:9)
    
    groupNames    = {'Pole', 'Lick', 'Reward', 'Pole + Lick', 'Pole + Reward', 'Lick + Reward', 'All factors', 'Other'};
    
%     plotStructMat  = zeros(4*8, 2);
%     
%     for nExample   = 1:8
%         for nPlot  = 1:5
%             plotStructMat(mod(nExample - 1, 4)*8+nPlot, ceil(nExample/4)) = nExample * 10 + nPlot;
%         end
%         plotStructMat(mod(nExample - 1, 4)*8+(6:7), ceil(nExample/4)) = nExample * 10 + 5;
%     end

    plotStructMat  = zeros(4*3, 2);
    
    for nExample   = 1:8
        for nPlot  = 1:3
            plotStructMat(mod(nExample - 1, 4)*3+nPlot, ceil(nExample/4)) = nExample * 10 + nPlot;
        end
    end


    plotStructMat = plotStructMat';
    
    unit_index    = zeros(4, minNumTrial);
    color_index   = [0.7 0.0 0.0; % CL
                     0.0 0.0 0.7; % CR
                     0.3 0.0 0.0; % EL
                     0.0 0.0 0.3]; % ER
                 
    whichNeuron   = [   3; % P
                        1; % L
                        3; % R
                        1; % PL
                        1; % PR
                        1; % LR
                        1; % PLR
                        1]; % other
    
    for nExample  = 1:8
        neuronIndex          = find(unitGroup == nExample & numTrialMat, whichNeuron(nExample));
        neuronIndex          = neuronIndex(whichNeuron(nExample));
        load([SpikingDataDir nDataSetList.cellinfo(neuronIndex).fileName], 'neuron_single_units', 'task_cue_time')
        currNeuronSpike      = neuron_single_units{nDataSetList.cellinfo(neuronIndex).nUnit};
        
        unit_index(1, :)     = nDataSet(neuronIndex).unit_yes_trial_index(1:minNumTrial);
        unit_index(2, :)     = nDataSet(neuronIndex).unit_no_trial_index(1:minNumTrial);
        unit_index(3, :)     = nDataSet(neuronIndex).unit_yes_error_index(1:minNumTrial);
        unit_index(4, :)     = nDataSet(neuronIndex).unit_no_error_index(1:minNumTrial);
        
        subplot(size(plotStructMat, 2), size(plotStructMat, 1), find(plotStructMat == (nExample * 10 + 1)));
        for nPlot            = 1:4
%             subplot(size(plotStructMat, 2), size(plotStructMat, 1), find(plotStructMat == (nExample * 10 + 1)));
            hold on;
            for nTrial       = 1:minNumTrial
                spikeTimes   = currNeuronSpike{unit_index(nPlot, nTrial)} - task_cue_time(unit_index(nPlot, nTrial));
                spikeTrial   = ones(length(spikeTimes), 1) * (nTrial + (nPlot-1)*(minNumTrial+3));
                plot(spikeTimes, spikeTrial, '.', 'color', color_index(nPlot, :));
            end
        end
        
        gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0);
        hold off;
        ylim([1 (minNumTrial + 3) * 4])
        xlim([params.timeSeries(1) params.timeSeries(end)]);
        title(groupNames{nExample});
        axis off
        
        subplot(size(plotStructMat, 2), size(plotStructMat, 1), find(plotStructMat == (nExample * 10 + 2)));
        hold on;
        nUnitData        = nDataSet(neuronIndex).unit_yes_trial;
        smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
        plot(params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(1, :));
        nUnitData        = nDataSet(neuronIndex).unit_no_trial;
        smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
        plot(params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(2, :));
        nUnitData        = nDataSet(neuronIndex).unit_yes_error;
        smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
        plot(params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(3, :));
        nUnitData        = nDataSet(neuronIndex).unit_no_error;
        smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
        plot(params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(4, :));
        gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
        hold off;
        legend('CA', 'CP', 'EA', 'EP', 'Orientation', 'horizontal', 'Location', 'northoutside');
        legend('boxoff');
        ylabel('Firing rate (Hz)');
        xlabel('Time (s)');
        xlim([params.timeSeries(1) params.timeSeries(end)]);
%         ylim([0 30]);
        
        
        subplot(size(plotStructMat, 2), size(plotStructMat, 1), find(plotStructMat == (nExample * 10 + 3)));
        hold on;
        h = plot(params.timeSeries, squeeze(logPValue(neuronIndex, :, :))','-', 'linewid', 1.0);
        gridxy ([params.polein, params.poleout, 0],[thresLogP], 'Color','k','Linestyle','--','linewid', 1.0)
        hold off;
        legend(h, 'Pole', 'Lick', 'Reward', 'Orientation', 'horizontal', 'Location', 'northoutside');
        legend('boxoff');
        ylabel('-log(P)');
        xlabel('Time (s)');
        xlim([params.timeSeries(1) params.timeSeries(end)]);
        ylim([0 maxLogP]);
        
    end
        
end

