# compute df/f from raw F data
#
# Ziqiang Wei
#
# weiz@janelia.hhmi.org

import numpy as np
from scipy.interpolate import interp1d


skip_interval = 3.0  # ignore F after a spike (< interval_pre sec) while computing baseline


def get_baseline_corr_dff(fmean, spktimes, ca_time, interval=skip_interval):
    fmeanSkipSpk = np.copy(fmean)
    for n_spk in spktimes:
        skip_frame = np.logical_and(ca_time >= n_spk, ca_time <= n_spk + interval)
        fmeanSkipSpk[skip_frame] = np.nan
    frameNoSpk = np.logical_not(np.isnan(fmeanSkipSpk))
    #  baseline = interp1d(ca_time[frameNoSpk], fmean[frameNoSpk], kind='linear')
    #  'cubic' interp1d behaves different in scipy from that in matlab, probably the num is different by default
    #  f_baseline = baseline(ca_time)
    f_baseline = fmean[frameNoSpk].mean()
    dff = (fmean - f_baseline) / f_baseline
    return dff
