function psth         = getPSTH(nData, trialType)
    averagedBinData   = mean(nData.(trialType), 1);
    boxCarWindowLength= 200; % ms
    boxCarWindow      = ones(1,boxCarWindowLength)/(boxCarWindowLength/1000);
    psth              = conv(averagedBinData, boxCarWindow, 'same');
    psth              = psth(201:end-200);
end