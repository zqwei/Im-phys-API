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


function [nDataSet3D, nDataSet] = getSimultaneousDataSet_v1_0(newDataSet, minUnitsSession)

    h                   = waitbar(0,'Initializing data analysis...');
    waitbar(0, h,'Low firing rate units are filtered out...');
    
    sessionIndex        = [newDataSet(:).sessionIndex];
    [sessionVec, ~, IC] = unique(sessionIndex);     
    valid_session       = hist(IC,length(sessionVec)) >= minUnitsSession;
    sessionVec          = sessionVec(valid_session);
    numSession          = length(sessionVec);   

    nDataSet            = [];
    nDataSet3D          = [];
    for nSession        = 1:numSession        
        nSessionData    = newDataSet(sessionIndex == sessionVec(nSession));
        [tSpikeDataSet, nodata, ~, nSessionData] = findInterSect(nSessionData);
        if ~nodata
            nDataSet3D  = [nDataSet3D; tSpikeDataSet]; %#ok<AGROW>
            nDataSet    = [nDataSet; nSessionData]; %#ok<AGROW>
        end
        waitbar(nSession/numSession, h, sprintf('%d of %d files have been finished...',nSession, numSession));
    end
    close (h)
end