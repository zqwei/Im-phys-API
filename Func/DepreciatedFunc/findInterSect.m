%
% findInterSect.m
%
% This file sets up the basic direction information of the Ca++ imaging
% data and the spiking data
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 




function [nInterSectData, nodata, interSectVec, nSessionData]   = findInterSect(nSessionData)
    numUnit            = length(nSessionData);
    
    interSectVec.unit_yes_trial_index     = intersect(nSessionData(1).unit_yes_trial_index, ...
                                                  nSessionData(2).unit_yes_trial_index);

    interSectVec.unit_no_trial_index      = intersect(nSessionData(1).unit_no_trial_index, ...
                                                  nSessionData(2).unit_no_trial_index);
                                              
    for nInterSect     = 3:numUnit
        interSectVec.unit_yes_trial_index = intersect(interSectVec.unit_yes_trial_index, ...
                                                      nSessionData(nInterSect).unit_yes_trial_index);        
        interSectVec.unit_no_trial_index  = intersect(interSectVec.unit_no_trial_index, ...
                                                      nSessionData(nInterSect).unit_no_trial_index);         
    end
    
    nInterSectData.sessionIndex           = nSessionData(1).sessionIndex;
    nInterSectData.nUnit                  = [nSessionData(:).nUnit];    
    nInterSectData.unit_yes_trial_index   = interSectVec.unit_yes_trial_index;
    nInterSectData.unit_no_trial_index    = interSectVec.unit_no_trial_index;
    if isempty(nInterSectData.unit_yes_trial_index) || isempty(nInterSectData.unit_no_trial_index)
        nodata = true;
    else
        nodata = false;
    end
    
    for nUnit          = 1:numUnit
        nSessionData(nUnit).unit_yes_trial         = nSessionData(nUnit).unit_yes_trial(...
                                                     ismember(nSessionData(nUnit).unit_yes_trial_index, nInterSectData.unit_yes_trial_index), :);
        nSessionData(nUnit).unit_yes_trial_index   = nInterSectData.unit_yes_trial_index;
        nInterSectData.unit_yes_trial(nUnit, :, :) = nSessionData(nUnit).unit_yes_trial;
        
        nSessionData(nUnit).unit_no_trial          = nSessionData(nUnit).unit_no_trial(...
                                                     ismember(nSessionData(nUnit).unit_no_trial_index, nInterSectData.unit_no_trial_index), :);                                         
        nSessionData(nUnit).unit_no_trial_index    = nInterSectData.unit_no_trial_index;
        nInterSectData.unit_no_trial(nUnit, :, :)  = nSessionData(nUnit).unit_no_trial;
    end
    
end
