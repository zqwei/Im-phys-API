%
% getLogPValueSpikeTime.m
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

function unitGroup     = getLogPValueTscoreSpikeTime(nDataSet, params)

    numUnit            = length(nDataSet);
    unitGroup          = zeros(numUnit, 1);
    minLength          = 5;
    
    for nUnit          = 1:numUnit
        % find p_value_string
        pValueTime    = getLogPValueTscore(nDataSet(nUnit), params);
        [M, V]        = regexp(sprintf('%i',[0 diff(pValueTime)==0] & pValueTime~=0),'1+','match');
        if ~isempty(M)
            lengthString = cellfun('length',M);
            V         = V(lengthString >= minLength);
            if ~isempty(V)
                unitGroup(nUnit) = length(unique(pValueTime(V)));
            end
        end
    end
    
end


function pValueTime    = getLogPValueTscore(spikeDataSet, params)
    sigma                         = 0.15 / params.binsize; % 300 ms
    filterLength                  = 11;
    filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
    filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
    filterInUse                   = filterInUse / sum (filterInUse); 

    nUnitData        = spikeDataSet.unit_yes_trial;
    yesUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
    nUnitData        = spikeDataSet.unit_no_trial;
    noUnitData       = getGaussianPSTH (filterInUse, nUnitData, 2);
    
    pValueTime       = ttest2(yesUnitData, noUnitData);
    signValueTime    = sign(mean(yesUnitData) - mean(noUnitData));
    pValueTime       = pValueTime .* signValueTime;
    
end