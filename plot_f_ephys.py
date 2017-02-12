import numpy as np
import h5py
from get_baseline_corr_dff import get_baseline_corr_dff
import matplotlib.pylab as plt

CaIndicator = 'GP43'
Zoom = 'highzoom' #  'highzoom'
SessionId = 'data_141001_cell1_001' #  'data_141001_cell1_001'

hf = h5py.File(CaIndicator + '_' + Zoom + '_nwb/' + SessionId + '.nwb', 'r')

# Raw fluorescence from nwb file
FMeanNeuropil = np.array(hf['acquisition']['timeseries']['FMeanNeuropil']['data'])
FMeanROI = np.array(hf['acquisition']['timeseries']['FMeanROI']['data'])
CaTimeStamps = np.array(hf['acquisition']['timeseries']['FMeanROI']['timestamps'])
print CaTimeStamps.max()

# Raw spike times from nwb file
SpkTimes = np.array(hf['processing']['Ephys']['UnitTimes']['SpikeTimes_0']['times'])

FMean = FMeanROI - 0.7 * FMeanNeuropil
MaxF = FMean.max()
MaxTime = CaTimeStamps.max()

# plot F vs spike times
fig, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(CaTimeStamps, FMean, '-b', label='Raw fluorescence')
axarr[0].plot(SpkTimes, np.ones(len(SpkTimes)) * (MaxF + 1.), '+k', label='Spike times')
axarr[0].legend()
axarr[0].axis([0, MaxTime, 0, MaxF + 5.])
axarr[0].set_xlabel('Time (s)')
axarr[0].set_ylabel('F')
axarr[0].set_title(CaIndicator + '_' + Zoom + '_' + SessionId)

# compute df/f
dff = get_baseline_corr_dff(FMean, SpkTimes, CaTimeStamps)
MaxDFF = dff.max()

# plot df/f vs spike times
axarr[1].plot(CaTimeStamps, dff, '-b', label='dF/F')
axarr[1].plot(SpkTimes, np.ones(len(SpkTimes)) * (MaxDFF + .5), '+k', label='Spike times')
axarr[1].legend()
axarr[1].axis([0, MaxTime, -1., MaxDFF + 1.])
axarr[1].set_xlabel('Time (s)')
axarr[1].set_ylabel('dF/F')
axarr[1].set_title(CaIndicator + '_' + Zoom + '_' + SessionId)
plt.show()
