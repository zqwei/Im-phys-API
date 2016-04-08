% 
% plot the spike dataset from a list of files
% 
% version 1.0
%
% Comparison list
%
% Output:
% SpikeDataSet     --- yDim x 1 cells (yDims number of neurons) 
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 



function plotSimultaneousSpikeData(nDataSet, params, minRate, perMinRate, ROCThres, minUnitsSession)
    newDataSet             = filterOutLowFR(nDataSet, params, minRate, perMinRate, ROCThres);
    trialUnitRatio         = 2.5;    
    sessionIndex           = unique([newDataSet.sessionIndex]);
    
    numSession             = 0;
    
    for nSession           = 1:length(sessionIndex)
        sessionUnitIndex   = find([newDataSet.sessionIndex] == sessionIndex(nSession));
        if length(sessionUnitIndex) >= minUnitsSession
            trialUnitMat   = zeros(32, 600);
%             figure; hold on;
            for nUnit          = 1:length(sessionUnitIndex)
%                 plot(newDataSet(sessionUnitIndex(nUnit)).unit_yes_trial_index,...
%                     newDataSet(sessionUnitIndex(nUnit)).nUnit*ones(length...
%                     (newDataSet(sessionUnitIndex(nUnit)).unit_yes_trial_index),1), ...
%                      'ok');
%                 plot(newDataSet(sessionUnitIndex(nUnit)).unit_no_trial_index,...
%                     newDataSet(sessionUnitIndex(nUnit)).nUnit*ones(length...
%                     (newDataSet(sessionUnitIndex(nUnit)).unit_no_trial_index),1), ...
%                      'or');
                trialUnitMat(newDataSet(sessionUnitIndex(nUnit)).nUnit, newDataSet(sessionUnitIndex(nUnit)).unit_yes_trial_index) = 1;
                trialUnitMat(newDataSet(sessionUnitIndex(nUnit)).nUnit, newDataSet(sessionUnitIndex(nUnit)).unit_no_trial_index) = -1;
            end
            
            trialIndex             = sum(trialUnitMat~=0, 1) >= minUnitsSession;
            unitIndex              = sum(trialUnitMat(:, trialIndex)==1, 2)> 20 & sum(trialUnitMat(:, trialIndex)==-1, 2)> 20;
            
            if sum(unitIndex) >= minUnitsSession;
%                 figure;                
%                 pcolor (trialUnitMat(unitIndex, trialIndex));
%                 yy = hist(sum(trialUnitMat~=0, 1), 1:15);
%                 plot(2:15, yy(2:15),'ok')
                trialIndex     = sum(trialUnitMat~=0, 1) == sum(unitIndex);
%                 trialUnitMat(unitIndex, trialIndex) = trialUnitMat(unitIndex, trialIndex) * 2;
%                 imagesc(trialUnitMat)
%                 title(sessionIndex(nSession))
                if sum(trialIndex & sum(trialUnitMat,1)<0) < sum(unitIndex) * trialUnitRatio ...
                    && sum(trialIndex & sum(trialUnitMat,1)>0) < sum(unitIndex) * trialUnitRatio
%                     disp([sessionIndex(nSession), sum(unitIndex), ...
%                          sum(trialIndex & sum(trialUnitMat,1)<0),  ...
%                          sum(trialIndex & sum(trialUnitMat,1)>0)])
%                     figure;
%                     imagesc(trialUnitMat)
%                     title(sessionIndex(nSession))
                    % kick one unit out
                    newTrialUnitMat     = trialUnitMat(unitIndex, :);
                    kickUnit            = 0;
                    numTrialAfterKick   = 0;
                    trialIndexAfterKick = false(1, size(newTrialUnitMat, 2));
                    for nKickUnit   = 1:size(newTrialUnitMat)
                        subTrialUnitMat = newTrialUnitMat;
                        subTrialUnitMat(nKickUnit, :) = [];
                        trialIndex     = sum(subTrialUnitMat~=0, 1) == size(subTrialUnitMat, 1);
                        
                        if numTrialAfterKick < sum(trialIndex)
                            numTrialAfterKick   = sum(trialIndex);
                            kickUnit            = nKickUnit;
                            trialIndexAfterKick = trialIndex;
                        end
                    end
                    
                    if kickUnit > 0
                        subTrialUnitMat = newTrialUnitMat;
                        subTrialUnitMat(kickUnit, :) = [];

                        if sum(trialIndexAfterKick & sum(subTrialUnitMat,1)<0) > size(subTrialUnitMat, 1) * trialUnitRatio ...
                            && sum(trialIndexAfterKick & sum(subTrialUnitMat,1)>0) > size(subTrialUnitMat, 1) * trialUnitRatio
                                disp([sessionIndex(nSession), sum(unitIndex)-1, ...
                                    sum(trialIndexAfterKick & sum(subTrialUnitMat,1)<0),  ...
                                    sum(trialIndexAfterKick & sum(subTrialUnitMat,1)>0)])
                                numSession = numSession + 1;
                        else
                            disp([sessionIndex(nSession), sum(unitIndex)-1, ...
                                    sum(trialIndexAfterKick & sum(subTrialUnitMat,1)<0),  ...
                                    sum(trialIndexAfterKick & sum(subTrialUnitMat,1)>0), -1])
                            figure;
                            imagesc(trialUnitMat)
                            title(sessionIndex(nSession))
                        end
                    else
                        disp([sessionIndex(nSession), sum(unitIndex), ...
                            sum(trialIndex & sum(trialUnitMat,1)<0),  ...
                            sum(trialIndex & sum(trialUnitMat,1)>0), -2]);
                        figure;
                        imagesc(trialUnitMat)
                        title(sessionIndex(nSession))
                    end
                else 
                    disp([sessionIndex(nSession), sum(unitIndex), ...
                        sum(trialIndex & sum(trialUnitMat,1)<0),  ...
                        sum(trialIndex & sum(trialUnitMat,1)>0)]);
                    numSession = numSession + 1;
                end
            end
            
        end
    end
    disp(numSession)
end

    

function HighFRDataSet = filterOutLowFR(nDataSet, params, minRate, perMinRate, ROCThres)
    
    HighFRDataSetIndex = arrayfun(@(x) mean(mean([x.unit_yes_trial; x.unit_no_trial])>minRate)>perMinRate, ...
                         nDataSet, 'Uniformoutput', false);
    HighFRDataSetIndex = cell2mat(HighFRDataSetIndex);
    
    ROCIndex           = ROCPop(nDataSet, params);
    ROCIndex           = ROCIndex > ROCThres;
    ROCDataIndex       = sum(ROCIndex(:,2:end),2)>0;
    
    HighFRDataSet      = nDataSet(HighFRDataSetIndex & ROCDataIndex);
    
end


