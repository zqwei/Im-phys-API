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

function coeffs       = coeffPCA(nSessionData)
    
    T                 = size(nSessionData, 3);

    coeffs            = arrayfun(@(tIndex) coeff1stPCA(nSessionData(:,:,tIndex)), 1:T, 'UniformOutput', false);
                                       
    coeffs            = cell2mat(coeffs);
    
end


function coeff        = coeff1stPCA(nSessionData)
    coeff             = pca(nSessionData,'numComponents',1);
    coeff             = coeff/norm(coeff);%coeff(:,1)./norm(coeff(:,1));
end