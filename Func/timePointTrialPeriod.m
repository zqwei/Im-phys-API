%
% timePointTrialPeriod.m
%
%
% ----------------------------
% Output:
% timePoints  --- [1 Polein PoleOut Cue End]
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function timePoints  = timePointTrialPeriod(polein, poleout, timeSeries)
    
    timePoints       = [1, sum(timeSeries<=polein), sum(timeSeries<=poleout)...
                        sum(timeSeries<=0), length(timeSeries)];