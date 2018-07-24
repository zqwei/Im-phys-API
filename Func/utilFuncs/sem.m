function y = sem(x)
    
    if size(x, 2) > 1
        y     = nanstd(x, [], 2);
        % len_y = sum(~isnan(y));
        len_y = length(x, 2);
        y     = y / sqrt(len_y);
    else
        y = nanstd(x, [], 1);
        % len_y = sum(~isnan(x));
        len_y = size(x, 1);
        y     = y / sqrt(len_y);
    end

end