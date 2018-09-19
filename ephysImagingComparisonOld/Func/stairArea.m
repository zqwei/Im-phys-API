%
% stairArea.m
%
% Plot stairwise area plot
% ----------------------------
% Output:
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function stairArea(x, y)

    x = [x;x];
    y = [y;y];
    area(x([2:end end]),y(1:end))

end