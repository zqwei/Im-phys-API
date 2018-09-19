function pValue = pValueTTest2 (x, y, filterInUse)

    if nargin        == 3
        x            = filter(filterInUse, 1, x, [], 2);
        y            = filter(filterInUse, 1, y, [], 2);
    end

    [~, pValue, ~]   = ttest2(x, y);

end