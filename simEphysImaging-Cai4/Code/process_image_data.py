# ------------------------------------------------------------------------------
# Function: create_time_series
# Returns: pointer to the created group
# Arguments: series_name, series_type, series_path, nwb_object,
#            var, t, description, source, comments,
#            hash_group , keyName, options
# Action:
# - creates a group <series_name>
#   of type <series_type>
#   at the path <spath>
#   relative to <nwb_group>
# - creates datasets "data", "num_samples" and "timestamps"
#   inside of the series group
# - populates the datasets from the arrays var and t
# - populates dataset "timestamps" with the data from input array t
#   using <keyName>
# - sets the attributes "description", "source" and "comments"
#   to the time series group using the input strings
#   <description>, <source> and <comments>
# - perform additional data handling using <h5_hash_group> and <keyName>
# NOTE: argument t has different meaning depending on the keyName
#       (it may be either start time or actual time)
def create_time_series(series_name, series_type, series_path, \
                       orig_h5, nwb_group, group_attrs, var, t, data_attrs, \
                       hash_group_pointer, keyName, options):
    # initialize timeseries
    if options.verbose:
        print "    Creating group ", series_name, " type=", series_type, " path=", series_path
    if len(series_path) > 0:
        ts = nwb_group.make_group(series_type, series_name, path=series_path, attrs=group_attrs)
    else:
        ts = nwb_group.make_group(series_type, series_name, attrs=group_attrs)
    if series_name in ["raw_ephys", "fmean_roi", "fmean_neurophil", \
                       "event_detected","filtered_ephys"]:
        timestamps = t
        on_off = np.int_(np.zeros(len(t)))
        on_off += -1
        on_off[::2] *= -1
        data = on_off
        # SP data only
        if series_name == "pole_touch_protract":
            kappa_ma_path = "eventSeriesArrayHash/value/2/eventPropertiesHash/1/value/2/2"
            kappa_ma = orig_h5[kappa_ma_path].value
            ts.set_custom_dataset("kappa_max_abs_over_touch", kappa_ma)
        elif series_name == "pole_touch_retract":
            kappa_ma_path = "eventSeriesArrayHash/value/2/eventPropertiesHash/2/value/2/2"
            kappa_ma = orig_h5[kappa_ma_path].value
            ts.set_custom_dataset("kappa_max_abs_over_touch", kappa_ma)
    # 2) Data = [1]*len(timestamps)
    elif series_name in ["lick_left", "lick_right",\
                         "pole_in",   "pole_out",  \
                         "lick_time"]              \
       or (series_name == "auditory_cue" and keyName == "CueTime"):
        data = [1] * len(t)
        timestamps = t
    # 3) Data = value
    elif series_name in ["whisker_angle", "whisker_curve", \
                         "lick_trace", "aom_input_trace",\
                         "simple_optogentic_stimuli"] \
       or keyName in ["whiskerVars", "Ephys"]:

        data = var
        timestamps = t
        if series_name == "simple_optogentic_stimuli":
            ts.set_dataset("site", "site 1")
    else:
       sys.exit("Unknown key "     + keyName     + \
                " or series_name " + series_name + \
                " in create_time_series")
    data_attrs['keyName'] = keyName
    ts.set_dataset("data", data, attrs=data_attrs)
    ts.set_dataset("timestamps", timestamps)
    ts.set_dataset("num_samples", len(timestamps))
    return ts

# ------------------------------------------------------------------------------
def save_2p_frames(external_file, starting_frame, timestamps, fname, stack_t):
    starting_frame.append(len(timestamps))
    timestamps.extend(stack_t)
    external_file.append(fname.encode('utf8'))

# ------------------------------------------------------------------------------
# pull out all ROI pixel maps for a particular subarea and imaging plane
#   and store these in the segmentation module
def fetch_rois(orig_h5, master_shape, plane_map, seg_iface, area, plane, \
               options, num_planes=3):
    tsah = orig_h5["timeSeriesArrayHash"]
    # convert from file-specific area/plane mapping to
    #   inter-session naming convention
    #image_plane = "area%d_plane%d" % (area, plane)
    oname = "area%d_plane%d" % (area, plane)
    master1_shape = master_shape[oname]
    image_plane = plane_map[oname]           

    # first get the list of ROIs for this subarea and plane
    if options.handle_errors:
        try:
            ids = tsah["value"]["%d"%(area+1)]["imagingPlane"]["%d"%plane]["ids"]
        except:
            # Warning: only one imaging plane is available (instead of 3)
            ids = tsah["value"]["%d"%(area+1)]["imagingPlane"]["ids"]
    else:
        ids = tsah["value"]["%d"%(area+1)]["imagingPlane"]["%d"%plane]["ids"]

    roi_ids = ids["ids"].value
    lookup = tsah["value"]["%d"%(area+1)]["ids"]["ids"].value
    for i in range(len(roi_ids)):
        rid = roi_ids[i]
        if num_planes == 1:
            rois = tsah["descrHash"]["%d"%(area+1)]["value"]["1"]
        else:
            rois = tsah["descrHash"]["%d"%(area+1)]["value"]["1"]["%d"%plane]
        # make sure the ROI id is correct
        try:
            record = rois["rois"]["%s"%(1+i)]
            x = int(parse_h5_obj(record["id"])[0])
            assert x == int(rid)
        except:
            print "Missing ROI for area=", area, " plane=", plane, " id=", i
            continue
        pix = parse_h5_obj(record["indicesWithinImage"])[0]
        # pix = record["indicesWithinImage/indicesWithinImage"].value
        pixmap = []
        for j in range(len(pix)):
            v = pix[j]
            px = int(v  / master1_shape[1])
            py = int(v) % master1_shape[0]
            pixmap.append([py,px])
        weight = np.zeros(len(pixmap)) + 1.0
#       print "image_plane=", image_plane, " oname=", oname
        ut.add_roi_mask_pixels(seg_iface, image_plane, "%d"%x, "ROI %d"%x, pixmap, \
            weight, master1_shape[1], master1_shape[0])

# ------------------------------------------------------------------------------
def fetch_dff(orig_h5, dff_iface, seg_iface, plane_map, area, plane, \
              options, num_planes=3):
    area_grp = orig_h5["timeSeriesArrayHash/value"]["%d"%(area+1)]
    if options.handle_errors:
        try:
            plane_ids = area_grp["imagingPlane"]["%d"%plane]["ids/ids"].value
        except:
            # Warning: only one imaging plane is available (instead of 3)
            plane_ids = area_grp["imagingPlane"]["ids/ids"].value
    else:
        plane_ids = area_grp["imagingPlane"]["%d"%plane]["ids/ids"].value

    area_ids = area_grp["ids/ids"].value
    # store dff in matrix and supply that to time series
    dff_data = area_grp["valueMatrix/valueMatrix"].value
    # convert from file-specific area/plane mapping to
    #   inter-session naming convention
    oname = "area%d_plane%d" % (area, plane)
    image_plane = plane_map[oname]
    t = area_grp["time/time"].value * 0.001
    # create array of ROI names for each matrix row
    roi_names = []
    trial_ids = area_grp["trial/trial"].value
    # for each plane ID, find group idx. df/f is idx'd row in values
    for i in range(len(plane_ids)):
        roi_names.append("ROI%d"%plane_ids[i])
    dff_ts = dff_iface.make_group("<RoiResponseSeries>", image_plane)
    dff_ts.set_dataset("data", dff_data, attrs={"unit":"dF/F",
        "conversion": 1.0, "resolution":0.0})
    dff_ts.set_dataset("timestamps", t)
    dff_ts.set_dataset("roi_names", roi_names)
    #- dff_ts.set_value_as_link("segmentation_interface", seg_iface)
    dff_ts.make_group("segmentation_interface", seg_iface)
    trial_ids = area_grp["trial/trial"].value
    dff_ts.set_custom_dataset("trial_ids", trial_ids)


# ------------------------------------------------------------------------------
def process_ROIs_and_dFoverF(orig_h5, dff_iface, seg_iface, num_subareas, \
                             plane_map, master_shape, options):
    # pull out image segmentation data. do it by subarea and imaging plane,
    #   as that's how data is stored in the source file
    print "Reading dF/F"
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        for plane in range(num_planes):
            sys.stdout.write('_')
    print ""
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        for plane in range(num_planes):
            fetch_rois(orig_h5, master_shape, plane_map, seg_iface, subarea+1, \
                       plane+1, options, num_planes)
            fetch_dff( orig_h5, dff_iface, seg_iface, plane_map, subarea+1, \
                       plane+1, options, num_planes)
            sys.stdout.write('.')
            sys.stdout.flush()
    print ""
    print "Writing ROI and dF/F"



# ------------------------------------------------------------------------------
def process_image_data(orig_h5, nwb_object, plane_map, options):
    # process dff
    mod = nwb_object.make_group("<Module>", "ROIs")
    mod.set_custom_dataset("description", "Segmentation (pixel-lists) and dF/F (dffTSA) for all ROIs")
    dff_iface = mod.make_group("DfOverF")
    dff_iface.set_attr("source", "This module's ImageSegmentation interface")
    seg_iface = mod.make_group("ImageSegmentation", attrs={"source": "Simon's datafile"})
    process_ROIs_and_dFoverF(orig_h5, dff_iface, seg_iface, num_subareas, \
                             plane_map, master_shape, options)


    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        for plane in range(num_planes):
            sys.stdout.write('_')
    print ""
    for subarea in range(num_subareas):
        # fetch time array
        grp = orig_h5["timeSeriesArrayHash"]["value"]["%d"%(subarea+2)]
        t = 0.001 * grp["time"]["time"].value
        # now move into imaging plane group, to generate time series for
        #   2photon image stacks
        grp = grp["imagingPlane"]
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        for plane in range(num_planes):
            if options.handle_errors:
                try:
                    pgrp = grp["%d"%(plane+1)]
                except:
                    # Warning: only one imaging plane is available (instead of 3)
                    pgrp = grp
            else:
                pgrp = grp["%d"%(plane+1)]

            frame_idx = pgrp["sourceFileFrameIdx"]["sourceFileFrameIdx"].value
            lst = parse_h5_obj(pgrp["sourceFileList"])[0]
            cnt = 0
            srcfile = {}
            for k in range(len(lst)):
                srcfile[str(k+1)] = lst[k]
                cnt += 1
            filemap = {}
            lastfile =-1
            lastframe = 1
            stack_t = []
            nname = None
            fname = None
            zero = np.zeros(1)
            assert len(t) == len(frame_idx[0])
            # following arrays used to make external_file as an array, reducing number of image_series
            external_file = []
            starting_frame = []
            timestamps = []
            for i in range(len(frame_idx[0])):
                filenum = frame_idx[0][i]
                if lastfile < 0:
                    lastfile = filenum
                framenum = frame_idx[1][i]
                stack_t.append(t[i])
                # check for embedded NaNs
                if np.isnan(filenum):
                    continue
                if np.isnan(framenum):
                    continue
                # use fname as a flag. if it's not None then there's data
                #   to write
                if fname is None:
                    # convert from file-specific area/plane mapping to
                    #   inter-session naming convention
                    oname = "area%d_plane%d" % (subarea+1, plane+1)
                    nname = plane_map[oname]
                    name = "%s_%d" % (nname, filenum)
                    fname = srcfile["%d"%filenum]
                # make sure frames and file numbers are sequential
                if not (lastfile == filenum and framenum == lastframe+1) and \
                   not                          framenum == 1:
                    # Warning: framenum or filenum does not start from 1
                    continue
#               assert (lastfile == filenum and framenum == lastframe+1) or framenum == 1
                if lastfile != filenum:
                    if i>0:
                        if not np.isnan(frame_idx[0][i-1] ) and \
                           not np.isnan(frame_idx[1][i-1]):
                            if options.handle_errors:
                                try:
                                     save_2p_frames(external_file, starting_frame, \
                                                    timestamps, fname, stack_t)
                                except:
                                    print "Warning: unable to create_2p_ts for name=", name
                            else:
                                save_2p_frames(external_file, starting_frame, \
                                               timestamps, fname, stack_t)
                            stack_t = []
                            fname = None
                lastframe = framenum
                lastfile = filenum
            # make sure we write out the last entry
            if fname is not None:
                if options.handle_errors:
                    try:
                        save_2p_frames(external_file, starting_frame, \
                                       timestamps, fname, stack_t)
                    except:
                        print "Warning: unable to create_2p_ts for name=", name
                else:
                    save_2p_frames(external_file, starting_frame, \
                                   timestamps, fname, stack_t)
            create_2p_tsa(nwb_object, nname, external_file, starting_frame, \
                          timestamps, name)
            sys.stdout.write('.')
            sys.stdout.flush()