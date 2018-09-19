
CaIndicator = 'GP43';
Zoom        = 'highzoom';
SessionId   = 'data_141001_cell1_001';



hinfo = hdf5info([CaIndicator '_' Zoom '_nwb/' SessionId '.nwb']);

fMeanNeuropil  = hdf5read(hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Datasets(1));
caTime         = hdf5read(hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Datasets(3));
fMeanROI       = hdf5read(hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(2).Datasets(1));
rawEphys       = hdf5read(hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(3).Datasets(1));
ephysTime      = hdf5read(hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(3).Datasets(4));
filteredEphys  = hdf5read(hinfo.GroupHierarchy.Groups(5).Groups.Groups(1).Groups.Datasets(1));
spkTimeEphys   = hdf5read(hinfo.GroupHierarchy.Groups(5).Groups.Groups(2).Groups.Datasets(2));

figure; % ephys
plot(ephysTime, filteredEphys)
hold on
plot(spkTimeEphys, ones(size(spkTimeEphys))*0.6, 'ro')

figure; % imaging
rawF           = fMeanROI - 0.7 * fMeanNeuropil;
plot(caTime, rawF)
hold on
plot(spkTimeEphys, ones(size(spkTimeEphys))*(max(rawF) + 5), 'ro')