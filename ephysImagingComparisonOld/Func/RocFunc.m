%
% plotROCExampleUnits.m
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

function [tp, fp]    = RocFunc(nRocTData, nRocOData, nBins)

    tPoints          = linspace(min(nRocTData), max(nRocTData), nBins);
    nYes             = sum(nRocOData==1);
    nNo              = sum(nRocOData==0);
    tp               = arrayfun(@(x) sum(x>nRocTData & nRocOData==1)/nYes, tPoints);
    fp               = arrayfun(@(x) sum(x>nRocTData & nRocOData==0)/nNo, tPoints);
end