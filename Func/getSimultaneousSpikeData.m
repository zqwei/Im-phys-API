% 
% obtain the spike dataset from a list of files
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



function [nDataSet3D, nDataSet] = getSimultaneousSpikeData(nDataSet, params, minRate, perMinRate, ROCThres, minUnitsSession)
    newDataSet             = filterOutLowFR(nDataSet, params, minRate, perMinRate, ROCThres);
    [nDataSet3D, nDataSet] = getSimultaneousDataSet(newDataSet, minUnitsSession);    
% % 
% %     h                   = waitbar(0,'Initializing data analysis...');
%     
%     newDataSet          = filterOutLowFR(nDataSet, params, minRate, perMinRate, ROCThres);
%     
%     waitbar(0, h,'Low firing rate units are filtered out...');
%     
%     sessionIndex        = [newDataSet(:).sessionIndex];
%     [sessionVec, ~, IC] = unique(sessionIndex);     
%     valid_session       = hist(IC,length(sessionVec)) >= minUnitsSession;
%     sessionVec          = sessionVec(valid_session);
%     numSession          = length(sessionVec);   
% %     SpikeDataSet        = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
% %                                 'unit_yes_trial', 1, 'unit_no_trial', 1,...
% %                                 'unit_yes_trial_index', 1, 'unit_no_trial_index', 1),numSession, 1);
%     
%     
%     tSession         = 1;
%     for nSession     = 1:numSession        
%         nSessionData           = newDataSet(sessionIndex == sessionVec(nSession));
%         [tSpikeDataSet, nodata] = findInterSect(nSessionData);
%         if ~nodata
%             SpikeDataSet(tSession) = tSpikeDataSet; %#ok<AGROW>
%             tSession               = tSession + 1;
%         end
%         waitbar(nSession/numSession, h, sprintf('%d of %d files have been finished...',nSession, numSession));
%     end
%     
%     close (h)
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

% function ROCIndex       = ROCPop(nDataSet, params)
% 
%     timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);     
%     numPlots            = length(nDataSet);
%     ROCIndex            = zeros(numPlots, length(timePoints)-1);
%     nBins               = 10;
%     for nPeriods        = length(timePoints) -1 : -1 : 1
%         nPeriodData     = dataInPeriods(nDataSet, timePoints, nPeriods);     
%         for nPlot       = 1:numPlots
%             nRocTData        = [nPeriodData(nPlot).unit_yes_trial; nPeriodData(nPlot).unit_no_trial];
%             nRocOData        = [ones(size(nPeriodData(nPlot).unit_yes_trial)); zeros(size(nPeriodData(nPlot).unit_no_trial))];
%             [tp, fp]         = RocFunc(nRocTData, nRocOData, nBins);
%             areaInt          = intXY([tp, 1], [fp, 1]);
%             ROCIndex(nPlot, nPeriods) = max(areaInt, 1 - areaInt);
%         end
%     end
%     
% end
% 
% 
% function areaInt       = intXY(vec_x, vec_y)    
%     areaInt            = sum((vec_y(1:end-1)+vec_y(2:end)).*(vec_x(2:end)-vec_x(1:end-1))/2);    
% end


% function [nInterSectData, nodata, interSectVec]   = findInterSect(nSessionData)
%     numUnit            = length(nSessionData);
%     
%     interSectVec.unit_yes_trial_index     = intersect(nSessionData(1).unit_yes_trial_index, ...
%                                                   nSessionData(2).unit_yes_trial_index);
% 
%     interSectVec.unit_no_trial_index      = intersect(nSessionData(1).unit_no_trial_index, ...
%                                                   nSessionData(2).unit_no_trial_index);
%                                               
%     for nInterSect     = 3:numUnit
%         interSectVec.unit_yes_trial_index = intersect(interSectVec.unit_yes_trial_index, ...
%                                                       nSessionData(nInterSect).unit_yes_trial_index);        
%         interSectVec.unit_no_trial_index  = intersect(interSectVec.unit_no_trial_index, ...
%                                                       nSessionData(nInterSect).unit_no_trial_index);         
%     end
%     
%     nInterSectData.sessionIndex           = nSessionData(1).sessionIndex;
%     nInterSectData.nUnit                  = [nSessionData(:).nUnit];    
%     nInterSectData.unit_yes_trial_index   = interSectVec.unit_yes_trial_index;
%     nInterSectData.unit_no_trial_index    = interSectVec.unit_no_trial_index;
%     if isempty(nInterSectData.unit_yes_trial_index) || isempty(nInterSectData.unit_no_trial_index)
%         nodata = true;
%     else
%         nodata = false;
%     end
%     
%     for nUnit          = 1:numUnit
%         nInterSectData.unit_yes_trial(nUnit, :, :) = nSessionData(nUnit).unit_yes_trial(...
%                                                      ismember(nSessionData(nUnit).unit_yes_trial_index, nInterSectData.unit_yes_trial_index), :);
%         nInterSectData.unit_no_trial(nUnit, :, :)  = nSessionData(nUnit).unit_no_trial(...
%                                                      ismember(nSessionData(nUnit).unit_no_trial_index, nInterSectData.unit_no_trial_index), :);
%     end
%     
% end
