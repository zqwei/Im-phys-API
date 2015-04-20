function z = zScore (x, y)

    z = (mean(x) - mean(y))./(var(x)/size(x,1) + var(y)/size(y,1));

end