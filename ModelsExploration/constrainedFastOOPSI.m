function constrainedFastOOPSI(dff, para)
    normalized_dff               = (dff - min(dff))/(max(dff) - min(dff));
    SAMP                         = cont_ca_sampler(normalized_dff);
    plot_continuous_samples(SAMP, normalized_dff, para);
