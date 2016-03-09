import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import fast_oopsi

T = 2000
dt = 0.020
lam = 0.1
tau = 1.5
sigma = 0.2

# fix random seed
np.random.seed(42)

# signal generator
F, C, N = fast_oopsi.fcn_generate(T, dt=dt, lam=lam, tau=tau, sigma=sigma)

# create figure
fig = plt.figure(num=None, figsize=(5, 7), dpi=90, facecolor='w',
                 edgecolor='k')
# fig = plt.figure()

# create subplot
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)

# operate on axis, plot grund-truth
ax1.plot(F)
# ax1.plot(C,color='red',linewidth=2)
ax2.plot(F, linewidth=2)
ax3.plot(F, linewidth=2)
ax4.plot(N, linewidth=2)
ax1.set_title('Synthetic Calcium Fluorescence Signal')

# fast-oopsi,
d, Cz, P, F_est = fast_oopsi.fast(F, dt=dt, iter_max=6)
ax2.plot(F_est, color='red', linewidth=1.5)
ax2.set_title('Reconstruct Spikes by fast-oopsi')

# wiener filter,
d, Cw, F_est, F_est_nonneg = fast_oopsi.wiener(F, dt=dt, iter_max=100)
ax3.plot(F_est_nonneg, color='red', linewidth=1.5)
ax3.set_title('Reconstruct Spikes by wiener filter')

# descritize,
d, v = fast_oopsi.discretize(F, bins=[0.75])
ax4.plot(d, color='red', linewidth=1.5)
ax4.set_title('Reconstruct Spikes by discretized binning')

# tunning the plot and show !
ax1.get_xaxis().set_visible(False)
ax2.get_xaxis().set_visible(False)
ax3.get_xaxis().set_visible(False)
ax1.set_yticks([-1, 0, 1])
ax2.set_yticks([0, 1])
ax3.set_yticks([0, 1])
ax4.set_yticks([0, 1])
font = {'weight': 'bold', 'size': 9}
matplotlib.rc('font', **font)
plt.show()
