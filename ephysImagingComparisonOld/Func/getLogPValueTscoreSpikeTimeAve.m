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

function unitGroup     = getLogPValueTscoreSpikeTimeAve(nDataSet, params, ffactor)

    numUnit            = length(nDataSet);
    unitGroup          = zeros(numUnit, 1);
    minLength          = 5;
    
    for nUnit          = 1:numUnit
        % find p_value_string
        pValueTime    = getLogPValueTscore(nDataSet(nUnit), params, ffactor);
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


function pValueTime  = getLogPValueTscore(spikeDataSet, params, ffactor)
    sigma            = 0.15 / params.binsize; % 300 ms
    filterLength     = 11;
    filterStep       = linspace(-filterLength / 2, filterLength / 2, filterLength);
    filterInUse      = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
    filterInUse      = filterInUse / sum (filterInUse); 
    
    

    nUnitData        = spikeDataSet.unit_yes_trial;
    yesUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
    nUnitData        = spikeDataSet.unit_no_trial;
    noUnitData       = getGaussianPSTH (filterInUse, nUnitData, 2);
    mean_yesUnitData = mean(yesUnitData, 1);
    mean_noUnitData  = mean(noUnitData, 1);
    yesSize          = length(spikeDataSet.unit_yes_trial_index);
    noSize           = length(spikeDataSet.unit_no_trial_index);
    numT             = size(nUnitData, 2);
    
    mean_rate        = min([mean_yesUnitData, mean_noUnitData]);
    baseline_rate    = 0;
%     ffactor          = 0.1;% 0.7;
    
    % std_yesUnitData  = std(yesUnitData) * ffactor;
    std_yesUnitData  = sqrt(mean_yesUnitData - mean_rate + baseline_rate) * ffactor;
    yesUnitData      = mean_yesUnitData + ones(size(yesUnitData,1), 1) * std_yesUnitData .* randn(yesSize, numT);
    
    
    % std_noUnitData   = std(noUnitData)* ffactor;
    std_noUnitData   = sqrt(mean_noUnitData - mean_rate + baseline_rate) * ffactor;
    noUnitData       = mean_noUnitData + ones(size(noUnitData, 1), 1) * std_noUnitData .* randn(noSize, numT);
    
    
    pValueTime       = ttest2(yesUnitData, noUnitData);
    signValueTime    = sign(mean(yesUnitData) - mean(noUnitData));
    pValueTime       = pValueTime .* signValueTime;
    
end