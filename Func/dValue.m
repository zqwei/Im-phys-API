function d = dValue (x, y)

    d = (mean(x) - mean(y))./max(abs(mean(x) - mean(y)));

end