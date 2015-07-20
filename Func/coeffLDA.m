%
% coeffLDA.m
%
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

function coeffs       = coeffLDA(nSessionData, totTargets)
    
    T                 = size(nSessionData, 3);

    coeffs            = arrayfun(@(tIndex) coeffClassify(...
                                           squeeze(nSessionData(:,:,tIndex)), ...
                                           totTargets), 1:T, 'UniformOutput', false);
                                       
    coeffs            = cell2mat(coeffs);
    
end


function coeff        = coeffClassify(varargin)

    obj               = fitcdiscr(varargin{:},'discrimType','pseudoLinear');
    coeff             = obj.Coeffs(1,2).Linear;
    coeff             = coeff./norm(coeff);
end