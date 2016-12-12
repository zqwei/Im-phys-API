function sigSelective = getSelective(nData, periodIndex)
    sigSelective      = false(1,length(periodIndex) - 2);
    for nPeriod       = 2:length(periodIndex)-1
        yesData       = mean(nData.unit_yes_trial(:,periodIndex(nPeriod):periodIndex(nPeriod+1)),2);
        noData        = mean(nData.unit_no_trial(:,periodIndex(nPeriod):periodIndex(nPeriod+1)),2);
        h             = ttest2(yesData, noData);
        if isnan(h); h = false; end
        sigSelective(nPeriod-1) = h;
    end
end