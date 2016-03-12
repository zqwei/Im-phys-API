import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import fast_oopsi
import constrained_oopsi
from scipy.io import loadmat, savemat


def KS_dat_fast_oopsi(nCell=1):
	totCell = loadmat('DataListCells.mat')
	dff = totCell['totCell'][nCell-1]['dff']
	dff = dff[0].astype('float64')[:,0]
	spk = totCell['totCell'][nCell-1]['spk']
	spk = spk[0]
	caTime = totCell['totCell'][nCell-1]['CaTime']
	caTime = caTime[0]
	dt = caTime[1] - caTime[0]
	# fast-oopsi,
    fast.d, fast.C, fast.P, fast.F_est = fast_oopsi.fast(dff, dt=dt, iter_max=100)
	# wiener filter,
    wiener.d, wiener.C, wiener.F_est, wiener.F_est_nonneg = fast_oopsi.wiener(dff, dt=dt, iter_max=100)
    # descritize,
    discr.d, discr.v = fast_oopsi.discretize(dff, bins=[0.75])
    # constrained-foopsi
    methods = ['cvxpy', 'spgl1', 'debug', 'cvx']
    cf.c, cf.bl, cf.c1, cf.g, cf.sn, cf.spikes = constrained_oopsi.constrained_foopsi(dff, p=2, noise_range=[.25, .5], methods=methods[0])
    return fast, wiener, discr, cf
