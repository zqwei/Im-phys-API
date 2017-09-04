%
% coeffLDA.m
%
%
% ----------------------------
% Output:
%
% version 1.0
%
% based on coeffLDA v 2.0
%
% only used for simultaneously recording datasets
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function [coeffs, errs] = coeffErrLDA(nSessionData, totTargets)
    
    T                 = size(nSessionData, 3);

    coeffs            = arrayfun(@(tIndex) coeffClassify(...
                                           squeeze(nSessionData(:,:,tIndex)), ...
                                           totTargets), 1:T, 'UniformOutput', false);
                                       
    coeffs            = cell2mat(coeffs);
    
end


function coeff        = coeffClassify(varargin)
    
    
    % choice of parameters follows 'min-min' rul in Guo's paper
    Mdl                    = fitcdiscr(varargin{:}, 'SaveMemory', 'on', ...
                            'FillCoeffs','off');
    [err, gamma, delta, ~] = cvshrink(Mdl, 'gamma', Mdl.MinGamma:0.015:0.5, ...
                            'delta', 0:0.03:1, 'Verbose', 0);
    % min-min error rule
    minerr                 = min(min(err));
    [p, q]                 = find(err == minerr);
    idx                    = sub2ind(size(delta), p, q); % Convert from subscripts to linear indices    
    gammaDeltaSet          = [gamma(p) delta(idx)];
    gammaDeltaSet          = gammaDeltaSet(gammaDeltaSet(:, 2)...
                                == max(gammaDeltaSet(:, 2)),:);
    gammaDeltaSet          = gammaDeltaSet(gammaDeltaSet(:, 1)...
                                == min(gammaDeltaSet(:, 1)),:);
    obj                    = fitcdiscr(varargin{:},'discrimType', 'Linear',...
                            'Delta', gammaDeltaSet(2), 'Gamma', gammaDeltaSet(1)); %pseudoLinear
                        
    coeff             = obj.Coeffs(1,2).Linear;
    coeff             = coeff./norm(coeff);
end