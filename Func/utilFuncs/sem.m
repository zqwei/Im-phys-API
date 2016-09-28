function y = sem(x)
    
    if size(x, 2) > 1
        y = std(x) / sqrt(size(x,2));
    else
        y = std(x) / sqrt(size(x,1));
    end

end