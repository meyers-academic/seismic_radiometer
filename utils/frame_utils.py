#! /usr/bin/env python
# Utilities for preparing data to be written to GWF frames using
# the Fr library.

#################################################################################
# Takes a dictionary of dataChannels which contain data for the entire day
# (data_dict) and makes a smaller dictionary containing data corresponding to the
# specified time period (frame_start_time to frame_start_time + frame_length).
# Also defines a DQ channel for the frame.
# INPUT:
#  data_dict        - dictionary of dictionaries. Entries correspond to station
#                     dictionaries. Each station dictionary contains dictionaries
#                     as well, where each of those dictionaries corresponds to
#                     data and metadata for a particular channel.
#  frame_start_time - GPS start time of the frame (all data samples should be
#                     at a time T >= frame_start_time.
#  frame_length     - duration of frame (seconds).
#  empty_data_val   - value corresponding to missing data.
#
# OUTPUT:
#  framedict_list - list of dictionaries where each entry corresponds to a
#                   particular station and channel. Each dictionary in the list
#                   contains data and metadata for writing frames. 
def compileFrameDict(data_dict,frame_start_time,frame_length,empty_data_val):

    # Imports.
    import numpy as np
    import sys, math

    # Calculate frame end (all samples should be LESS than this time).
    frame_end_time = frame_start_time + frame_length

    # List of dictionaries for all stations and channels.
    # Each entry will be a dictionary for a particular station/channel combination.
    framedict_list = []

    for sta_i in data_dict.keys():
        # List of dictionaries for all channels for a particular station.
        chandict_list = []
        for chan_j in data_dict[sta_i].keys():

            # Set up data object.
            dataObj = data_dict[sta_i][chan_j]

            # Get indices for data to add to channel dictionary.
            # Start index is determined based on difference between
            # frame start time and data object's start time.
            t_diff = frame_start_time - dataObj.startTime
            if (t_diff <= 0):
                if (t_diff >= dataObj.dt):
                    raise ValueError('Data start time is more than 1 dt less than ' + \
                                     'frame_start_time...something went wrong')
                else:
                    start_idx = 0
            else:
                start_idx = int(math.ceil(t_diff/dataObj.dt))

            # End index is based on frame length.
            end_idx = start_idx + (int(frame_length/dataObj.dt)-1)
            end_time = end_idx*dataObj.dt + dataObj.startTime
            # For channels with dt = 20, there can be either 6 or 7 samples within
            # a 128 second range, depending on when it starts; this block is designed to
            # account for this.
            if (end_time < (frame_end_time - dataObj.dt)):
                end_idx += 1
                end_time += dataObj.dt
            
            # Check start/end times to be sure they are within a sample time
            # and fall in the correct range.
            start_time = start_idx*dataObj.dt + dataObj.startTime
            if ((start_time-frame_start_time) > dataObj.dt or start_time < frame_start_time):
                raise ValueError('Data start time does not fall in the correct range.')
            if ((frame_end_time-end_time) > dataObj.dt or end_time >= frame_end_time):
                raise ValueError('Data end time does not fall in the correct range.')

            # Build dictionary for this channel.
            chandict = {}
            chandict['name'] = dataObj.station + ':' + dataObj.channel # Channel name.
            chandict['data'] = dataObj.data[start_idx:(end_idx+1)] # numpy array of channel data.
            chandict['start'] = frame_start_time # FRAME start time.
            # startX is time of first sample in channel RELATIVE to frame start time.
            chandict['startX'] = dataObj.startTime + start_idx*dataObj.dt - frame_start_time
            chandict['dx'] = dataObj.dt # sampling time.
            chandict['x_unit'] = 's' # X unit (sampling times)
            chandict['y_unit'] = 'counts' # Y unit (data units)
            chandict['kind'] = 'ADC' # not sure, but this is what Shivaraj used.
            chandict['type'] = 1 # 1 = time-series data.

            # Append chandict to list of channel dictionaries.
            chandict_list.append(chandict)

        # Append DQ channel.
        DQdict = getDQchannel(sta_i,chandict_list,frame_start_time,frame_length,empty_data_val)
        chandict_list.append(DQdict)

        # Add to full list of dictionaries for all station/channel combinations.
        framedict_list.extend(chandict_list)

    return framedict_list



#################################################################################
# Returns string which specifies the frame filename based on frame root directory,
# frame directory prefix and frame prefix, GPS start time, and duration.
# Path to frames: (frame_root)/(framedir_prefix)-XXXXX/(frame_prefix)-(GPS_start)-
#                 (dur).gwf, where XXXXX is the first 5 digits of the frame's
#                 GPS start time.
# INPUT:
#  frame_root      - path to root directory where frames will be saved.
#  framedir_prefix - prefix for directories where frames will be stored.
#  frame_prefix    - prefix for naming frames.
#  GPS_start       - GPS start time of the frame.
#  duration        - frame duration (seconds).
#
# OUTPUT:
#  frame_name      - full name of frame.
#  framedir        - path where frame will be stored.
def getFrameName(frame_root,framedir_prefix,frame_prefix,GPS_start,duration):
    # Imports
    import os

    # Assumes convention of (framedir_prefix)-XXXXX/frame_prefix-(GPS_start)-(dur).gwf
    # where XXXXX is the first 5 digits of the GPS start time.
    framedir_num = str(GPS_start)[:5]
    framedir = frame_root + '/' + framedir_prefix + '-' + str(framedir_num)

    # frame name
    frame_name = framedir + '/' + frame_prefix + '-' + \
        str(GPS_start) + '-' + str(int(duration)) + '.gwf'

    return frame_name, framedir

#################################################################################
# Checks various data quality channels to define an overall DQ channel
# for the frame.  The result is a single number which is between 0 and
# 31 depending on which of 5 flags are found to be true.
# INPUT:
#  sta_name       - name of the station for which we calculating the DQ channel.
#  chandict       - dictionary containing data objects for all channels for
#                   this station.
#  frame_start    - GPS start time of the frame for which we are calculating
#                   the DQ channel.
#  frame_length   - duration of the frame (seconds).
#  empty_data_val - value corresponding to missing data.
#
# OUTPUT:
#  DQdict - dictionary for this station's  DQ channel (single number indicating
#           the quality of the station's data in this frame). The dictionary
#           contains all of the necessary metadata for being written to a frame.
def getDQchannel(sta_name,chandict,frame_start,frame_length,empty_data_val):

    # imports
    import numpy as np
    from itertools import groupby, count
    from homestake_meta import getChannels

    # Useful channels for DQ flag:
    # LCQ - clock quality
    # VM1 - mass position for channel 1
    # VM2 - mass position for channel 2
    # VM3 - mass position for channel 3
    #
    # Data is encoded as one number per station per frame.
    # This value is a set of binary values:
    # 2^0: clock quality dips below 100% at some point
    # 2^1: mass position for channel 1 goes outside the range [-40,40]
    # 2^2: mass position for channel 2 goes outside the range [-40,40]
    # 2^3: mass position for channel 3 goes outside the range [-40,40]
    # 2^4: some or all data is missing
    # If the value of the DQ channel is 0, that means everything is awesome!

    # Threshold on consecutive data elements equal to missing_data_val
    # for data to be considered missing.
    missing_thresh = 5

    # Clock threshold: require 100% clock quality.
    clk_thresh = 100
    
    # Mass positions - should be between [-40,40] counts (or [-4,4] volts).
    mp_low = -40
    mp_high = 40

    # Values for DQ channel.
    dq_chans = {}
    dq_chans['LCQ'] = 1
    dq_chans['VM1'] = 2
    dq_chans['VM2'] = 4
    dq_chans['VM3'] = 8
    dq_chans['AMD'] = 16 # AMD = "Any Missing Data"

    # Set up DQ dictionary.
    DQdict = {}
    DQdict['name'] = sta_name + ':ODQ' # ODQ = "Overall Data Quality"
    DQdict['data'] = 0
    DQdict['start'] = frame_start
    DQdict['startX'] = frame_start
    DQdict['dx'] = frame_length
    DQdict['x_unit'] = 's'
    DQdict['y_unit'] = 'counts'
    DQdict['kind'] = 'ADC' # not sure, but this is what Shivaraj used
    DQdict['type'] = 1 # 1 = time-series


    # Check data quality. --------------
    # Check clock quality.
    clk_data = [chan['data'] for chan in chandict if chan['name'].find('LCQ') >= 0][0]
    if (min(clk_data) < clk_thresh):
        DQdict['data'] += dq_chans['LCQ']

    # Check mass position 1.
    mp1_data = [chan['data'] for chan in chandict if chan['name'].find('VM1') >= 0][0]
    if (min(mp1_data) < mp_low or max(mp1_data) > mp_high):
        DQdict['data'] += dq_chans['VM1']

    # Check mass position 2.
    mp2_data = [chan['data'] for chan in chandict if chan['name'].find('VM2') >= 0][0]
    if (min(mp2_data) < mp_low or max(mp2_data) > mp_high):
        DQdict['data'] += dq_chans['VM2']

    # Check mass position 3.
    mp3_data = [chan['data'] for chan in chandict if chan['name'].find('VM3') >= 0][0]
    if (min(mp3_data) < mp_low or max(mp3_data) > mp_high):
        DQdict['data'] += dq_chans['VM3']

    # Only use this next part for data stations (i.e., everything except MAST).
    if (sta_name.find('MAST') == -1):
        # Get data channels.
        data_chan = getChannels('data')
        # Check each one for missing data.
        missing_data = False
        for chan_i in data_chan:
            temp_data = [chan['data'] for chan in chandict if chan['name'].find(chan_i) >= 0][0]
            
            # Find groups of missing_data_val.
            # We assume that elements corresponding to missing_data_val
            # can occur randomly but probably not consecutively.
            missing_idx = np.nonzero(temp_data == empty_data_val)[0].tolist()
            
            if (len(missing_idx) == 0):
                consec_missed = 0
            else:
                groups = groupby(missing_idx, key=lambda item, c=count():item-next(c))
                sizes = [list(g) for k, g in groups]
                consec_missed = max([len(x) for x in sizes])
                
                if (consec_missed > missing_thresh):
                    DQdict['data'] += dq_chans['AMD']
                    break

    # Make data a numpy array so it matches the format of other channels.
    DQdict['data'] = np.asarray(DQdict['data'])

    return DQdict

# EOF
