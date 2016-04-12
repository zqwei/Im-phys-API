
function noiseCorrelationEpochs
    addpath('../Func');
    setDir;
    load ([TempDatDir 'DataListShuffle.mat']);
    for nData   = [1 3 4]
        load([TempDatDir DataSetList(nData).name '.mat'])
        [yesTrialSignalCorr, yesTrialNoiseCorr, noTrialSignalCorr, noTrialNoiseCorr, yes_no_pair]...
            = computeSignalCorrAndNoiseCorrEpochs(nDataSet, DataSetList(nData).params);
        
        figure;
        titles = {'Sample', 'Delay', 'Response'};
        for nPlot = 1:3
            subplot(4, 3, nPlot)
            [f1, xi1] = hist(yesTrialSignalCorr(:, nPlot), 1000);
            [f2, xi2] = hist(noTrialSignalCorr(:, nPlot), 1000);
            plot(xi1, f1/sum(f1), '-b', xi2, f2/sum(f2), '-r', 'linewid', 1.0)
            ylabel('Frac.')
            xlabel('Signal corr.')
            box off
            title(titles{nPlot})
            
            subplot(4, 3, nPlot+3)
            [f1, xi1] = hist(yesTrialNoiseCorr(:, nPlot), 1000);
            [f2, xi2] = hist(noTrialNoiseCorr(:, nPlot), 1000);
            plot(xi1, f1/sum(f1), '-b', xi2, f2/sum(f2), '-r', 'linewid', 1.0)
            ylabel('Frac.')
            xlabel('Noise corr.')
            box off
            
            subplot(4, 3, nPlot+3*2)
            plot(yesTrialSignalCorr(yes_no_pair(:, 1), nPlot), noTrialSignalCorr(yes_no_pair(:, 2), nPlot), '.k', [-1 1], [-1 1], '--r', 'linewid', 1.0)
            ylabel('Signal corr. (ipsi.)')
            xlabel('Signal corr. (contra.)')
            box off

            
            subplot(4, 3, nPlot+3*3)
            plot(yesTrialNoiseCorr(yes_no_pair(:, 1), nPlot), noTrialNoiseCorr(yes_no_pair(:, 2), nPlot), '.k', [-1 1], [-1 1], '--r', 'linewid', 1.0)
            ylabel('Noise corr. (ipsi.)')
            xlabel('Noise corr. (contra.)')
            box off

            
        end
        
        setPrint(8*3, 6*4, [PlotDir 'SingleUnitsNoiseCorrelation\' DataSetList(nData).name])
        setPrint(8*3, 6*4, [PlotDir 'SingleUnitsNoiseCorrelation\' DataSetList(nData).name],'png')
        
    end
end

function [yesTrialSignalCorr, yesTrialNoiseCorr, noTrialSignalCorr, noTrialNoiseCorr, yes_no_pair]...
            = computeSignalCorrAndNoiseCorrEpochs(nDataSet, params)
    
    numUnit            = length(nDataSet);
    yesTrialSignalCorr = nan(numUnit, 3);
    yesTrialNoiseCorr  = nan(numUnit, 3);
    noTrialSignalCorr  = nan(numUnit, 3);
    noTrialNoiseCorr   = nan(numUnit, 3);
    yes_no_pair        = nan(numUnit, 2);
    timePoints         = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);   
    
    nPairYes           = 0;
    nPairNo            = 0;
    nPairYesNo         = 0;
    allSessionInd      = [nDataSet.sessionIndex];
    sessionInd         = unique(allSessionInd);
    
    for nSession       = 1:length(sessionInd)
        nSessInd       = sessionInd(nSession);
        if sum(allSessionInd==nSessInd) > 1 % session with pairs
            unitInd    = find(allSessionInd==nSessInd);
            for nUnit  = 1:length(unitInd)
                for mUnit = nUnit+1:length(unitInd)
                    % yes trial
                    % find pairs
                    nUnitIdx      = unitInd(nUnit);
                    mUnitIdx      = unitInd(mUnit);
                    nUnitTrialInd = nDataSet(nUnitIdx).unit_yes_trial_index;
                    mUnitTrialInd = nDataSet(mUnitIdx).unit_yes_trial_index;
                    [C1, ai, bi]   = intersect(nUnitTrialInd, mUnitTrialInd);
                    if length(C1) > 10
                        nPairYes  = nPairYes + 1;
                        yesTrialN = nDataSet(nUnitIdx).unit_yes_trial(ai, :);
                        yesTrialM = nDataSet(mUnitIdx).unit_yes_trial(bi, :);
                        for nPeriod = 1:3
                            yesTrialNPeriod = yesTrialN(:,timePoints(nPeriod+1):timePoints(nPeriod+2));
                            yesTrialMPeriod = yesTrialM(:,timePoints(nPeriod+1):timePoints(nPeriod+2));
                            yesTrialSignalCorr(nPairYes, nPeriod) = corr(mean(yesTrialNPeriod, 1)', mean(yesTrialMPeriod, 1)');
                            yesTrialNoiseCorr(nPairYes, nPeriod)  = corr(mean(yesTrialNPeriod, 2), mean(yesTrialMPeriod, 2));
                        end
                    end
                    
                    % no trial
                    % find pairs
                    nUnitIdx      = unitInd(nUnit);
                    mUnitIdx      = unitInd(mUnit);
                    nUnitTrialInd = nDataSet(nUnitIdx).unit_no_trial_index;
                    mUnitTrialInd = nDataSet(mUnitIdx).unit_no_trial_index;
                    [C2, ai, bi]   = intersect(nUnitTrialInd, mUnitTrialInd);
                    if length(C2) > 10
                        nPairNo  = nPairNo + 1;
                        noTrialN = nDataSet(nUnitIdx).unit_no_trial(ai, :);
                        noTrialM = nDataSet(mUnitIdx).unit_no_trial(bi, :);
                        for nPeriod = 1:3
                            noTrialNPeriod = noTrialN(:,timePoints(nPeriod+1):timePoints(nPeriod+2));
                            noTrialMPeriod = noTrialM(:,timePoints(nPeriod+1):timePoints(nPeriod+2));
                            noTrialSignalCorr(nPairNo, nPeriod) = corr(mean(noTrialNPeriod, 1)', mean(noTrialMPeriod, 1)');
                            noTrialNoiseCorr(nPairNo, nPeriod)  = corr(mean(noTrialNPeriod, 2), mean(noTrialMPeriod, 2));
                        end
                    end
                    
                    if length(C1)>10 && length(C2)>10
                        nPairYesNo         = nPairYesNo + 1;
                        yes_no_pair(nPairYesNo, :) = [nPairYes, nPairNo];
                    end
                end
            end
        end
    end
    
end