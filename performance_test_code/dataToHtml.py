from glob import glob
import pandas as pd
from bokeh.palettes import Category20
from bokeh.plotting import figure, output_file, show
colors= Category20.get(20)
names = ['dF/F', 'S2C Linear', 'S2C Hill', 'S2C Sigmoid', 'C2S NWF', 'C2S FOOPSI', 'C2S FRI', 'C2S AR1', 'C2S AR2', 'C2S AR3', 'C2S MCMC', 'C2S Peel']


for nFolder in glob('MV1*'):
    spkFiles = glob(nFolder + '/*_spk.csv')
    dffFiles = glob(nFolder + '/*_dff.csv')
    for nFile in range(len(spkFiles)):
        fileName = spkFiles[nFile][:-8]#[len(nFolder)+1:]
        outputFileName = fileName + '.html'
        output_file(outputFileName)
        spkdb = pd.read_csv(spkFiles[nFile])
        spkdb['spk'] = spkdb['spk'] * -1.0
        dffdb = pd.read_csv(dffFiles[nFile])
        p = figure(plot_width=1400, plot_height=500, x_axis_label="Time(sec)", y_axis_label="dF/F")
        p.circle(spkdb['spkTime'], spkdb['spk'], line_width=2, color=colors[0], legend='Spikes')
        for name, color in zip(names, range(1, len(names)+1)):
            p.line(dffdb['time'], dffdb[dffdb.columns[color]], line_width=2, color=colors[color], legend=name)
        p.legend.location = "top_left"
        p.legend.click_policy="hide"
        p.legend.orientation = "horizontal"
        show(p)
