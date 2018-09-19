function spkCount     = getSpkCount(nData, spkCountTime, spkCountDur, trialType)
    averagedBinData   = mean(nData.(trialType), 1);
    spkCount          = sum(averagedBinData(spkCountTime))/spkCountDur;
end