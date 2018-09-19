%
% getGaussianPSTH.m
% 
%
% Spiking dataset
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

function smoothedUnitData = getGaussianPSTH (filterInUse, nUnitData, fDim)

    if fDim ~=2
        error('Fileter dimension is not supported');
    end
    
    zeroPadLength = floor(length(filterInUse)/2);
    numData       = size(nUnitData, 1);
    zeroPaddedData = [zeros(numData, zeroPadLength), nUnitData, zeros(numData, zeroPadLength)];
    weightMat     = [zeros(numData, zeroPadLength), ones(size(nUnitData)), zeros(numData, zeroPadLength)];
    
    smoothedUnitData = filter(filterInUse, 1, zeroPaddedData, [], fDim);
    weightMat      = filter(filterInUse, 1, weightMat, [], fDim);
    
%     smoothedUnitData = smoothedUnitData(:, zeroPadLength+1:end-zeroPadLength)./ ...
%                         weightMat(:, zeroPadLength+1:end-zeroPadLength);

    smoothedUnitData = smoothedUnitData(:, 1+2*zeroPadLength:end)./ ...
                        weightMat(:, 1+2*zeroPadLength:end);%weightMat(:, zeroPadLength+1:end-zeroPadLength);

end