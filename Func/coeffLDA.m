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

function coeffs       = coeffLDA(nSessionData, trainingTargets, testTargets)
    
    T                 = size(nSessionData, 3);

    coeffs            = arrayfun(@(tIndex) coeffClassify(...
                                           squeeze(nSessionData(1:length(testTargets),:,tIndex)), ...
                                           squeeze(nSessionData(length(testTargets)+1:end,:,tIndex)), ...
                                           trainingTargets), 1:T, 'UniformOutput', false);
                                       
    coeffs            = cell2mat(coeffs);
    
end


function coeff        = coeffClassify(varargin)
    
    [~,~,~,~,coeff]   = classify(varargin{:});
    coeff             = coeff(1,2).linear;
    coeff             = coeff./norm(coeff);
    
end