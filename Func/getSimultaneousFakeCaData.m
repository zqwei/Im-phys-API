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



function [nDataSet3D, nDataSet] = getSimultaneousFakeCaData(nCaDataSet, nSpikeDataSet, params, minRate, perMinRate, ROCThres, minUnitsSession)
    newDataSet             = filterOutLowFR(nCaDataSet, nSpikeDataSet, params, minRate, perMinRate, ROCThres);
    [nDataSet3D, nDataSet] = getSimultaneousDataSet(newDataSet, minUnitsSession);    
%     
% 
%     h                   = waitbar(0,'Initializing data analysis...');
%     
%     newDataSet          = filterOutLowFR(nCaDataSet, nSpikeDataSet, params, minRate, perMinRate, ROCThres);
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

    

function HighFRDataSet = filterOutLowFR(nCaDataSet, nSpikeDataSet, params, minRate, perMinRate, ROCThres)
    
    HighFRDataSetIndex = arrayfun(@(x) mean(mean([x.unit_yes_trial; x.unit_no_trial])>minRate)>perMinRate, ...
                         nSpikeDataSet, 'Uniformoutput', false);
    HighFRDataSetIndex = cell2mat(HighFRDataSetIndex);
    
    ROCIndex           = ROCPop(nSpikeDataSet, params);
    ROCIndex           = ROCIndex > ROCThres;
    ROCDataIndex       = sum(ROCIndex(:,2:end),2)>0;
    
    HighFRDataSet      = nCaDataSet(HighFRDataSetIndex & ROCDataIndex);
    
end