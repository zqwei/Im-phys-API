function d = dPrime (x, y)

    d = (mean(x) - mean(y))./(mean(x) + mean(y));

end