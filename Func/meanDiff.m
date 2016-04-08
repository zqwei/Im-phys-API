function meanDiffValue = meanDiff (x, y, filterInUse)

    if nargin        == 3
        x            = filter(filterInUse, 1, x, [], 2);
        y            = filter(filterInUse, 1, y, [], 2);
    end

    meanDiffValue    = mean(x) - mean(y);

end