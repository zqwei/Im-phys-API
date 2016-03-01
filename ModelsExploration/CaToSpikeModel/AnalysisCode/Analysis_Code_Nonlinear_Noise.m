tCell                        = 40;
para.t_frame                 = totCell(tCell).CaTime;
para.t_ephys                 = totCell(tCell).ephysTime;
para.fneuropil               = totCell(tCell).fNeuropil;
para.fmean                   = totCell(tCell).fROI;
para.filt                    = totCell(tCell).filteredEphys;
para.peak                    = totCell(tCell).detectedSpikes;
dff                          = get_baseline(para);
% dff                          = para.fmean - 0.7 * para.fneuropil;

binSize    = 10000;
spikeCouts = sum(reshape(para.peak, binSize, []), 1);
sum(spikeCouts == 1)
spikeTime = para.t_ephys((1:length(spikeCouts))*binSize);
spike_To_Plot = find(spikeCouts == 1);
totCaTime     = [];
figure;
hold on
for iPlot = 1:length(spike_To_Plot)
    nSpikeTime = min(para.t_ephys(para.t_ephys >= spikeTime(spike_To_Plot(iPlot)) & para.peak == 1));
    if ~isempty(nSpikeTime)
        nCaTime    = find(para.t_frame >= nSpikeTime, 1);
        totCaTime     = [totCaTime; nCaTime];
        plot((1:100)/60, dff(nCaTime:nCaTime+99),'k');
    end
end
ylabel('df/f')
xlabel('Time after spike (s)')
setPrint(8, 6, 'Fig5', 'pdf')

figure;
hold on;
plot(para.t_frame, dff);
plot(para.t_frame(totCaTime), 5, '+r')
hold off
ylabel('df/f')
xlabel('Time (s)')
setPrint(16, 6, 'Fig6', 'pdf')