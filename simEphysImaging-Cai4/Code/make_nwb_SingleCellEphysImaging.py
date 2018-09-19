#!~/anaconda/bin/python
#
# make_nwb.py
#
# Created by Ziqiang Wei on 2016-09-03.

import sys
import os

# import nwb
from nwb import nwb_file
from nwb import nwb_utils as ut

import h5py
import datetime
# import getpass
import numpy as np
# from sets import Set
import re
import optparse
import h5lib

os.environ['NWB_DATA'] = './'  # path of "surgery.txt", "data_collection.txt", "experiment_description.txt"


# ------------------------------------------------------------------------------
def make_nwb_command_line_parser(parser):
    parser.add_option("-D", "--debug", action="store_true", dest="debug", help="output debugging info", default=False)
    parser.add_option("-e", "--no_error_handling", action="store_false", dest="handle_errors", help="handle_errors", default=True)
    parser.add_option("-o", "--outfolder", dest="output_folder", help="output folder (default=same as input folder)", metavar="output_folder", default="")
    parser.add_option("-r", "--replace", action="store_true", dest="replace", help="if the output file already exists, replace it", default=False)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="increase the verbosity level of output", default=False)

    return parser


# ------------------------------------------------------------------------------
def check_entry(file_name, obj):
    try:
        return file_name[obj]
    except KeyError:
        print str(obj) + " does not exist"
        return []


# ------------------------------------------------------------------------------
# Extract all keys from input data and check if there are keys not in the pre-defined list
def check_keys(orig_h5, meta_h5, options):
    known_keys = ['indicator', 'timeSeriesArrayHash', 'timeUnitNames', 'trialTimeUnit']
    if options.verbose:
        print "known_keys=", sorted(known_keys)
    unknown_keys = []
    all_keys = h5lib.get_all_keys(orig_h5, meta_h5)
    if options.verbose:
        print "\nAll keys=", sorted(all_keys)
    for k in all_keys:
        if k not in known_keys and \
           not re.search("unit", k) and \
           not re.search("dffTSA", k):
            unknown_keys.append(k)
    return


# ------------------------------------------------------------------------------
def parse_h5_obj(obj, level=0, output=[], verbose=0):
    if level == 0:
        output = []
    try:
        if isinstance(obj, h5py.highlevel.Dataset):
            level = level + 1
            if obj.value.any():
                output.append(obj.value)
            else:
                full_name = obj.name.split("/")
                output.append(full_name[-1])
        elif isinstance(obj, h5py.highlevel.Group):
            level = level + 1
            if not obj.keys():
                output.append([])
            else:
                for key in obj.keys():
                    parse_h5_obj(obj[key], level, output, verbose)
        else:
            output.append([])
    except KeyError:
        print "Can't find" + str(obj)
        output.append([])
    return output


# ------------------------------------------------------------------------------
def set_metadata(group, keyname, value):
    print "    set_metadata: keyname=", keyname, " value=", value
    if keyname in ["extracellular", "intracellular"]:
        group.set_dataset("description", value)
    else:
        print "keyname=", keyname, " value=", value
        group.set_custom_dataset(keyname, value)


# ------------------------------------------------------------------------------
def set_metadata_from_file(group, keyname, file_name):
    group.set_dataset(keyname, ut.load_file(os.path.join(os.environ['NWB_DATA'], file_name)))


# ------------------------------------------------------------------------------
def process_metadata(nwb_object, input_h5, session_id, options):
    # input_h5 is from meta_h5 file
    if options.verbose:
        print "Processing metadata"
    general_group = nwb_object.make_group("general", abort=False)
    subject_group = general_group.make_group("subject", abort=False)
    # for (key, value) in zip(["age", "weight"], [0, 0]):
    #     set_metadata(subject_group, key, value)
    # session_id: 'data_'
    #           + datetime (5:10)
    #           + 'cell#' (12:16)
    #           + repetition (18:20)
    set_metadata(subject_group, "cell", session_id[12:16])
    # general_group.set_custom_dataset("reference_atlas", value)
    # for key in ["virus", "fiber", "photostim", "surgicalManipulation"]:
    #     general_group.set_custom_dataset(key, value)
    # set_metadata(subject_group, "subject_id", value) # aminal ID
    set_metadata(general_group, "session_id", session_id)
    genotype = 'Thy1-GCaMP6-WPRE GP5.17'
    set_metadata(subject_group, "genotype", genotype + "\n")
    # set_metadata(subject_group, "description", "animalStrain:" + animalStrain + "  animalSource:" + animalSource + "  dateOfBirth:" + subject "\n")
    set_metadata(general_group, "notes", "Simultaneous recordings of intracellular calcium levels and action potentials of individual neurons in the visual cortex during visual stimulations of drifting gratings.")
    set_metadata(general_group, "related_publications", "--")
    set_metadata(subject_group, "sex", "male")
    set_metadata(subject_group, "species", "mouse")
    set_metadata_from_file(general_group, "surgery", "surgery.txt")
    set_metadata_from_file(general_group, "data_collection", "data_collection.txt")
    set_metadata_from_file(general_group, "experiment_description", "experiment_description.txt")
    set_metadata(general_group, "experimenter", "Bei-Jung Lin")
    set_metadata(general_group, "institution", "Janelia Research Campus, HHMI")
    set_metadata(general_group, "lab", "Svoboda's lab")


# ------------------------------------------------------------------------------
def process_image_data(orig_h5, nwb_object, options):
    # process fmean_roi data in acquisition
    vs = nwb_object.make_group("<TimeSeries>", "FMeanROI", path="/acquisition/timeseries")
    data = parse_h5_obj(orig_h5['timeSeriesArrayHash/value/1/valueMatrix'])[0]
    vs.set_attr("comments", "Acquired at 55Hz; 256 x 256 pixels")
    vs.set_attr("source", "Device 'imaging-acquisition'")
    vs.set_dataset("data", data, attrs={"unit": "a.u.", "conversion": 1.0, "resolution": float('nan')})
    timestamps = parse_h5_obj(orig_h5['timeSeriesArrayHash/value/1/time'])[0]
    vs.set_dataset("timestamps", timestamps, attrs={"unit": "Seconds"})
    vs.set_attr("description", "Time series to represent average fluorescence of ROI.")

    # process fmean_roi data in acquisition
    vs = nwb_object.make_group("<TimeSeries>", "FMeanNeuropil", path="/acquisition/timeseries")
    data = parse_h5_obj(orig_h5['timeSeriesArrayHash/value/2/valueMatrix'])[0]
    vs.set_attr("comments", "Acquired at 55Hz; 256 x 256 pixels")  # area of neuropil
    vs.set_attr("source", "Device 'imaging-acquisition'")
    vs.set_dataset("data", data, attrs={"unit": "a.u.", "conversion": 1.0, "resolution": float('nan')})
    timestamps = parse_h5_obj(orig_h5['timeSeriesArrayHash/value/2/time'])[0]
    vs.set_dataset("timestamps", timestamps, attrs={"unit": "Seconds"})
    vs.set_attr("description", "Time series to represent average fluorescence of surrounding neuropil.")


# ------------------------------------------------------------------------------
def process_extracellular_spike_time(orig_h5, nwb_object, options):
    # process raw ephys data in acquisition
    vs = nwb_object.make_group("<ElectricalSeries>", "RawEphys", path="/acquisition/timeseries")
    data = parse_h5_obj(orig_h5['timeSeriesArrayHash/value/3/valueMatrix'])[0]
    vs.set_attr("comments", " ")
    # vs.set_attr("comments", "Acquired at 19531.25 Hz")  # is the correct value?
    vs.set_attr("source", "Device 'ephys-acquisition'")
    vs.set_dataset("data", data, attrs={"unit": "Volts", "conversion": 0.001, "resolution": float('nan')})
    timestamps = parse_h5_obj(orig_h5['timeSeriesArrayHash/value/3/time'])[0]
    vs.set_dataset("timestamps", timestamps, attrs={"unit": "Seconds"})
    el_idx = [0]
    vs.set_dataset("electrode_idx", el_idx)
    vs.set_attr("description", "Time series to represent raw ephys recording.")

    # process filtered ephys
    mod = nwb_object.make_group("<Module>", "Ephys")
    fephys = mod.make_group("FilteredEphys")
    fephys.set_attr("source", "Data as reported in Bei-Jung's mat file")
    fephys = fephys.make_group("<ElectricalSeries>", 'unit_0')
    fephys.set_attr("comments", " ")
    # fephys.set_attr("comments", "Acquired at 19531.25 Hz")  # is the correct value?
    fephys.set_attr("description", 'Filtered ephys for each channel at 5kHz')
    fephys.set_attr("source", "Data as reported in Bei-Jung's mat file")
    data = parse_h5_obj(orig_h5['timeSeriesArrayHash/value/4/valueMatrix'])[0]
    timestamps = parse_h5_obj(orig_h5['timeSeriesArrayHash/value/4/time'])[0]
    fephys.set_dataset("data", data, attrs={"resolution": float('nan'), "unit": "Volts", "conversion": 0.001})
    fephys.set_dataset("timestamps", timestamps, attrs={"unit": "Seconds"})
    el_idx = [0]
    fephys.set_dataset("electrode_idx", el_idx)

    # process spike events
    fephys = mod.make_group("UnitTimes")
    fephys.set_attr("source", "Data as reported in Bei-Jung's mat file")
    # fephys.set_attr("description", "Times for detected spike events")
    fephys = fephys.make_group("<unit_N>", 'SpikeTimes_0')
    fephys.set_dataset("unit_description", "Times for detected spike events; digitized at 10kHz")
    fephys.set_dataset("source", "Data as reported in Bei-Jung's mat file")
    # fephys.set_attr("source_electricalseries", nwb_object['/acquisition/timeseries/RawEphys'])
    # fephys.set_attr("source_electricalseries_path", "/acquisition/timeseries/RawEphys")
    # fephys.set_attr("source_idx", [0])
    data = parse_h5_obj(orig_h5['timeSeriesArrayHash/value/5/valueMatrix'])[0]
    timestamps = parse_h5_obj(orig_h5['timeSeriesArrayHash/value/5/time'])[0]
    fephys.set_dataset("times", timestamps[data == 1])
    # fephys.set_dataset("data", data)
    # fephys.set_dataset("timestamps", timestamps, attrs={"unit": "Seconds", "conversion": 0.001, "resolution": float('nan')})


# ------------------------------------------------------------------------------
def produce_nwb(data_path, metadata_path, output_nwb, options):
    orig_h5 = h5py.File(data_path, "r")
    meta_h5 = ""
    if len(metadata_path) > 0:
        meta_h5 = h5py.File(metadata_path, "r")
    check_keys(orig_h5, meta_h5, options)
    print "Input:", data_path, " ", metadata_path
    print "Output:", output_nwb, "\n"
    if options.replace and os.path.exists(output_nwb):
        os.remove(output_nwb)
    # the keys follow description of data structure
    vargs = {}
    session_id = os.path.basename(output_nwb)[0:-4]
    # session_id: 'data_'
    #           + datetime (5:10)
    #           + 'cell#' (12:16)
    #           + repetition (18:20)
    #
    # experiment date
    dt = datetime.datetime.strptime("20" + session_id[5:11], "%Y%m%d")
    print(session_id[5:10])
    print(dt.strftime("%a %b %d %Y"))
    vargs["start_time"] = dt.strftime("%a %b %d %Y")
    vargs["file_name"] = output_nwb
    vargs["identifier"] = ut.create_identifier(session_id)
    # a description txt file for experiment
    vargs["description"] = 'Date : ' + dt.strftime("%a %b %d %Y") + '\n' + \
                    'Cell ID : ' + session_id[12:17] + '\n' + \
                    'Session ID: ' + session_id[18:21] + '\n' + \
                    'Zoom: High zoomed'
    # vargs["description"] = os.path.join(os.environ['NWB_DATA'],"experiment_description.txt")
    vargs["mode"] = "w"
    # print "vargs=", vargs
    nwb_object = nwb_file.open(**vargs)
    print "Processing metadata ..."
    # No "metaDataHash" in this set of data, process_metaData is using independent metaDataHash file, c.f. Nuo...
    process_metadata(nwb_object, meta_h5, session_id, options)
    # processing of imaging data
    print "Processing imaging data ..."
    process_image_data(orig_h5, nwb_object, options)
    # processing of ephys
    print "Processing extracellular ephys data ..."
    # process_ephys_electrode_map(nwb_object, meta_h5, options)
    process_extracellular_spike_time(orig_h5, nwb_object, options)
    # for future release of S2C model
    print "Processing analysis of S2C model ..."
    print "----skip"
    # collect_analysis_information(orig_h5, nwb_object, options)
    if options.verbose:
        print "Closing file"
    nwb_object.close()

# ------------------------------------------------------------------------------


if __name__ == "__main__":
    usage = "Usage: \n%prog data_h5 [meta_data_h5] [options (-h to list)]"
    parser = optparse.OptionParser(usage=usage)
    parser = make_nwb_command_line_parser(parser)
    (options, args) = parser.parse_args()
    if options.verbose:
        print "len(args)=", len(args)
    if len(args) in [1, 2]:
        data_path = args[0]
        if not re.search(".h5", data_path):
            sys.exit("\nScript make_mwb.py accepts as input .h5 data file")
        metadata_path = ""
        if len(args) == 2:
            metadata_path = args[1]
            if not re.search(".h5", metadata_path):
                sys.exit("\nScript make_mwb.py accepts as input .h5 metadata file")
        data_basename = os.path.basename(data_path)
        if len(options.output_folder) == 0:
            options.output_folder = os.path.dirname(data_path)
        output_path = os.path.join(options.output_folder, data_basename.split(".")[0] + ".nwb")
        produce_nwb(data_path, metadata_path, output_path, options)
    else:
        parser.print_usage()
        sys.exit(2)
