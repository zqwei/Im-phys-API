%
% explainedPCA.m
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

function normlizedData        = normalizationDim(nSessionData, xDim)

    totDim                    = ndims(nSessionData);
    reorderDim                = 1:totDim;
    reorderDim(xDim)          = [];
    reorderDim                = [reorderDim, xDim];    
    xLength                   = size(nSessionData, xDim);    
    xData                     = reshape(permute(nSessionData, reorderDim), [], xLength);
    meanData                  = mean(xData);
    stdData                   = std(xData);
    
    if xDim                   == 1
        normlizedData         = bsxfun(@minus, nSessionData, meanData');
        normlizedData         = bsxfun(@rdivide, normlizedData, stdData');
    else
        normlizedData         = bsxfun(@minus, nSessionData, meanData);
        normlizedData         = bsxfun(@rdivide, normlizedData, stdData);        
    end
    
    
    
end