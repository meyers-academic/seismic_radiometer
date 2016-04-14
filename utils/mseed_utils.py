#! /usr/bin/env python
# Miscellaneous functions for handling data, including
# loading miniSEED files and handling the data.

#################################################################################
# Used to calculate number of samples needed to span a time period time_in with
# sampling time dt. Useful due to annoying rounding issues in Python.  If
# the fractional number of samples is larger than 1/100 of sample, we round up,
# otherwise we round down.
# INPUT:
#  time_in - time period in seconds.
#  dt      - sampling time in seconds.
#
# OUTPUT:
#  nSamples - number of samples.
def getSamples(time_in, dt):
    nSamples = time_in/float(dt)
    delta = nSamples - int(nSamples)
    if (delta >= 1/100.0):
        nSamples = int(nSamples) + 1
    else:
        nSamples = int(nSamples)
    
    return nSamples

#################################################################################
# Get sampling time from channel name.
# INPUT:
#  channel_name - three-letter channel name (HHZ, LCQ, etc.).
#
# OUTPUT:
#  dt - sampling time (seconds).
def getDeltaT(channel_name):
    # Check channel name.  delta_t is defined as
    # H**: 0.01 s, L**: 1 s, V**: 10 s, Q**: 20 s
    fs_dict = {'H':0.01, 'L':1, 'V':10, 'Q':20}
    return fs_dict[channel_name[0]]

#################################################################################
# Get name of miniSEED file based on a hard-coded naming convention.
# INPUT:
#  db_root - the root path to the miniSEED database.
#  net     - network code (X6 for this experiment)
#  sta     - station name
#  chan    - channel code
#  UTCdate - UTCDateTime object of date of interest.
#
# OUTPUT:
#  fname_str - db_root/YYYY/DDD/station.net..chan.YYYY.DDD
#    (here YYYY is the 4 digit year and DDD is the 3 digit day of the year)
def getMSEEDfname(db_root,net,sta,chan,UTCdate):
    # Set up filename.
    fname_str = db_root + '/' + UTCdate.strftime('%Y/%j') + '/' + sta + '.' + \
        net + '..' + chan + '.' + UTCdate.strftime('%Y.%j')

    return fname_str


# Takes in obspy data structure (containing multiple blocks) and returns GPS times
# corresponding to the first and last data samples which span the day.
# INPUT:
#  data_array - obspy array of data blocks.
#  dt         - sampling time.
#
# OUTPUT:
#  t_first - GPS time of first sample.
#  t_last  - GPS time of last sample.
#
# Example: data at 0.1 Hz with sample times 12, 22, 32, ..., 86392. t_first will be 2,
# since there should have been a sample at that time.
def getDataTimeRange(data_array, dt, day_GPSstart, day_length):

    # Imports.
    import math
    from timing_utils import utc2gps

    # Loop over all blocks and get UTC times of first and last samples.
    utc_first = data_array[0].stats.starttime
    utc_last = data_array[0].stats.endtime
    for i in range(1, len(data_array)):
        if (data_array[i].stats.starttime < utc_first):
            utc_first = data_array[i].stats.starttime
        if (data_array[i].stats.endtime > utc_last):
            utc_last = data_array[i].stats.endtime
    # Convert to GPS, subtract day_GPSstart, and round to
    # nearest sampling time.
    # t_first and t_last are relative to day_GPSstart to handle
    # rounding errors.
    nDigits = max(int(math.floor(math.log10(1/float(dt)))), 0)
    first_sample = round(utc2gps(utc_first)-day_GPSstart, nDigits)
    last_sample = round(utc2gps(utc_last)-day_GPSstart, nDigits)

    # Check if the file contains data in the correct time range.
    if (last_sample < 0 or first_sample > day_length):
        # If not, return.
        data_flag = False
        return (None, None, False)
    else:
        # Set data_flag to True.
        data_flag = True

    # We want to return the time of the first and last samples that contain data
    # during the requested interval. For cases where dt <= 1, t_first should be
    # 0 (i.e., 0 seconds RELATIVE to day_GPSstart).
    if (dt <= 1):
        t_first = 0
    else:
        # For cases with dt > 1, the samples might occur at times which don't
        # intersect with day_GPSstart. Example: 0.1 Hz data sampled at -2, 8,
        # 18, 28, etc. We handle these cases here.   
        # Adjust times as needed.
        if (first_sample >= dt):
            # If first_sample is at least dt (after day_GPSstart).
            t_first = first_sample % dt
        elif (first_sample < 0):
            # If first_sample is less than 0 (i.e. earlier than day_GPSstart).
            #t_first = dt - (-first_sample % dt)
            t_first = (first_sample % dt)
        else:
            # Otherwise, t_first is just first_sample.
            t_first = first_sample

    # Time of last sample.
    t_last = t_first + day_length - dt

    return (t_first, t_last, data_flag)

#################################################################################
# Returns start and end times of the day in GPS time and leap second variables.
# INPUT:
#  day_UTCstart - UTCdatetime object representing the start time of the day.
#
# OUTPUT:
#  day_GPSstart - start time of the day in GPS seconds.
#  day_length   - length of day in seconds (86400 normally, 86401 for leap seconds).
#  ls_today     - True if there was a leap second added on this day.
#  ls_yest      - True if there was a leap second added on the previous day.
#
# Note: because day_GPSend is actually the start time of the next day, all data
# samples for today should be at a time T such that day_GPSstart <= T < day_GPSend.
def getDayInfo(day_UTCstart):
    # Imports.
    import datetime
    from timing_utils import utc2gps, getLeapSeconds

    # Here day_GPSstart is the start time of the current day and
    # day_GPSend is the start time of the next day. So all samples
    # should be at some time t >= day_GPS_start and t < day_GPSend.
    day_GPSstart = utc2gps(day_UTCstart)
    day_GPSend = utc2gps(day_UTCstart + datetime.timedelta(days=1))
    day_length = round(day_GPSend-day_GPSstart)

    # Check if a leap second is implemented today.
    ls_today = False
    if (day_length == 86401):
        ls_today = True

    # Check if there was a leap second implemented yesterday.
    _, ls_prev = getLeapSeconds(day_UTCstart + datetime.timedelta(days=-1))
    ls_yest = False
    if (ls_prev > 0):
        ls_yest = True

    return (day_GPSstart, day_length, ls_today, ls_yest)

#################################################################################
# Gets GPS times of first and last samples in a block.
# INPUT:
#  stats        - stats from an obspy block.
#  dt           - sampling time for this channel.
#  day_GPSstart - GPS time of the beginning of the day.
#
# OUTPUT:
#  block_tfirst - time of first sample in block relative to day_GPSstart.
#  block_tlast  - time of last sample in block relative to day_GPSstart.
def getBlockTimes(stats,dt,day_GPSstart):
    # Imports.
    import math
    from timing_utils import utc2gps

    # Start/end times from data block rounded to nearest second or dt,
    # depending on sampling frequency.
    nDigits = max(int(math.floor(math.log10(1/float(dt)))), 0)
    block_tfirst = round(utc2gps(stats.starttime)-day_GPSstart,nDigits)
    block_tlast = block_tfirst + dt*(stats.npts-1)

    return (block_tfirst, block_tlast)

#################################################################################
# Adjusts block start/end times for a leap second.
# The Q330s seems to handle leap seconds by starting a new block which overlaps
# with a previous block by 1 second.  The break may come slightly before the leap
# second or slightly after it, we check for both cases. Note that there is a
# hard-coded value of 30 seconds - this is the time before or after the leap second
# in which the break may occur.  This is based on my experience with the 6/30/2015
# leap second and may need to be adjusted in the future.
# INPUT:
#  data           - array of obspy data blocks.
#  i              - index of current block.
#  block_tfirst   - time of first sample in the block relative to day_GPSstart.
#  block_tlast    - time of last sample in the block relative to day_GPSstart.
#  day_GPSstart   - GPS start time of the day.
#  day_length     - 86400 normally; 86401 for days with a leap second.
#  ls_today       - boolean which tells whether there was a leap second on this day.
#  ls_yest        - boolean which tells whether there was a leap second on the previous day.
#
# OUTPUT:
#  block_tfirst - (see inputs).
#  block_tlast  - (see inputs).
def adjustTimesForLeapSecond(data,i,dt,block_tfirst,block_tlast,day_GPSstart,day_length,ls_today,ls_yest):
    # The Q330/Antelope seems to implement a break (i.e. a new block)
    # near the time of the leap second (either before or after).
    import pdb
    #pdb.set_trace()
    break_window = 30 # seconds

    # Check if the leap second was today.
    if ((ls_today and (day_length - block_tfirst) < break_window and i > 0) or
        (ls_yest and block_tfirst < 0 and block_tfirst > -break_window and i > 0)):
        # For cases where the break was before the leap second.
        if (abs(data[i].stats.starttime - data[i-1].stats.endtime) < dt/10.0):     
            block_tfirst += 1; block_tlast += 1
    elif (ls_yest and block_tlast < break_window and i == 0 and len(data) > 1):
        # If the leap second was yesterday, we check for cases which
        # have the break implemented after the leap second.
        #if (data[i].stats.endtime == data[i+1].stats.starttime or
        #    data[i].stats.endtime == data[i+1].stats.starttime - dt + 1):
        if (abs(data[i].stats.endtime - data[i+1].stats.starttime) < dt/10.0 or
            abs(data[i].stats.endtime - (data[i+1].stats.starttime-dt+1)) < dt/10.0):
            block_tfirst -= 1; block_tlast -= 1
    


    return (block_tfirst, block_tlast)

#################################################################################
# Gets indices of block data which fall between day_GPSstart and day_GPSend
# Returns -1 for all outputs if no relevant data is found.
# INPUT:
#  block_tfirst - time of first sample in block relative to day_GPSstart.
#  block_tlast  - time of last sample in block relative to day_GPSstart.
#  stats        - obspy block information.
#  dt           - sampling time (seconds).
#  day_GPSstart - GPS start time of the day.
#  day_length   - day length in seconds (usually 86400, but 86401 for days with a leap second).
#
# OUTPUT:
#  c_idx1      - index corresponding to first sample from the block which will be
#                falls between day_GPSstart and day_GPSend
#  c_idx2      - same as c_idx2, but for last sample in the block in the correct range.
#  data_tfirst - time of sample at c_idx1, relative to day_GPSstart.
#  data_tlast  - time of sample at c_idx2, relative to day_GPSstart.
def getBlockIndices(block_tfirst, block_tlast, stats, dt, day_GPSstart, day_length):
    # Imports.
    import math

    # Check if block contains any relevant data.
    if (block_tfirst >= day_length or block_tlast < 0):
        return (None,None,None,None)

    # Get index and GPS time of first sample in the block which falls between
    # day_GPSstart and day_GPSend.
    c_idx1 = int(math.ceil((max(-block_tfirst,0)/float(dt))))
    data_tfirst = block_tfirst + c_idx1*dt

    # Get index and GPS time of last sample in the block which falls between
    # day_GPSstart and day_GPSend.
    if (block_tlast >= day_length):
        # Case where block_tlast >= day_length.
        c_idx2 = getSamples((day_length-dt)-block_tfirst,dt)
    else:
        # Case where 0 <= block_tlast < day_length.
        c_idx2 = int(stats.npts - 1)
    data_tlast = data_tfirst + dt*(c_idx2-c_idx1)

    return (c_idx1, c_idx2, data_tfirst, data_tlast)

#################################################################################
# Gets indices corresponding to the full data array where we will copy the
# current block data into. I.e., chanObj.data[idx1:(idx2+1)] = data[i].data[c_idx1:(c_idx2+1)].
# INPUT:
#  data_tfirst - time of first sample in block relative to day_GPSstart.
#  t_first     - time of first data sample for the day.
#  t_last      - time of last data sample for the day.
#                (all times relative to the start time of the day).
#  dt          - sampling time (seconds).
#  nSamp       - number of samples to be copied from the block to the full data array.
#
# OUTPUT:
#  idx1 - index of the full data array where we will copy the first relevant block sample.
#  idx2 - index of the full data array where we will copy the last relevant block sample.
def getDataIndices(data_tfirst, t_first, t_last, dt, nSamp):

    # Get indices where we want to copy the data to.
    # nSamp is calculated from block indices.
    idx1 = getSamples((data_tfirst - t_first),dt)
    idx2 = idx1 + nSamp

    # Catch for weird cases with dt = 10 or 20.
    t2 = data_tfirst + (idx2-idx1)*dt
    if (t2 > t_last and dt > 1):
        idx1 -= 1
        idx2 -= 1

    return (idx1, idx2)

#################################################################################
# Function used to load the next day's data (if there is a leap second added on
# the current day) and add the leap second data to the current day's data. This is
# necessary because of the way the leap seconds are implemented in the Q330s;
# the leap second data may be saved in the next day's data.
# INPUT:
#  chanObj      - object containing channel data and metadata.
#  fname        - filename corresponding to the current day.
#  day_UTCstart - current day's start time as a UTCDateTime object.
# OUTPUT:
#  (chanObj is returned by reference)
def addLeapSecondDataFromTomorrow(chanObj,fname,day_UTCstart):
    # Imports
    from obspy.core import read
    import datetime, re, pdb
    import numpy as np

    # Get start time of tomorrow.
    tom_UTCstart = day_UTCstart + datetime.timedelta(days=1)

    # Set up filename and load data.
    m = re.search(r'^(.*?)\/\d{4}.*\.(\w{2})\.\.',fname)
    db_root = m.group(1); net = m.group(2)
    tomorrow_fname = getMSEEDfname(db_root,net,chanObj.station,chanObj.channel,tom_UTCstart)
    data = read(tomorrow_fname)

    # Loop over blocks, get first second of data.
    for i in range(0,len(data)):
        if (data[i].stats.starttime.second < 1):
            # Indices of data to pull from the block.
            b_idx1 = 0
            b_idx2 = int(((tom_UTCstart + 1) - data[i].stats.starttime)/float(chanObj.dt)) - 1

            # Indices of new_data to put the data in.
            d_idx1 = int((data[i].stats.starttime - day_UTCstart)/chanObj.dt)
            d_idx2 = d_idx1 + (b_idx2 - b_idx1)
            
            # Add to chanObj.data.
            chanObj.data[d_idx1:(d_idx2+1)] = data[i].data[b_idx1:(b_idx2+1)]

# Get data from miniSEED file.  Data will be stored in chanObj.data.
def getDataFromMSEED2(chanObj,day_UTCstart,empty_data_val):
    # Imports
    import os, math, datetime, sys
    import numpy as np
    from obspy.core import read
    from timing_utils import utc2gps, getLeapSeconds

    # Set up shorthand variables.
    dt = chanObj.dt
    nSamples = chanObj.nSamples

    # Get GPS times of start and end of day.
    # Here day_GPSstart is the start time of the current day and
    # day_GPSend is the start time of the next day. So all samples
    # should be at some time t >= day_GPS_start and t < day_GPSend.
    day_GPSstart = utc2gps(day_UTCstart)
    day_GPSend = utc2gps(day_UTCstart + datetime.timedelta(days=1))

    # Check if a leap second is implemented today.
    ls_today = False
    if (day_GPSend - day_GPSstart == 86401):
        ls_today = True

    # Check if there was a leap second implemented the previous day.
    _, ls_prev = getLeapSeconds(day_UTCstart + datetime.timedelta(days=-1))
    ls_yest = False
    if (ls_prev > 0):
        ls_yest = True

    # Load data.
    data = read(chanObj.fname)
        
    # Set up time array
    # Loop over all blocks and get earliest time
    utc_first = data[0].stats.starttime
    utc_last = data[0].stats.endtime
    for i in range(1,len(data)):
        if (data[i].stats.starttime < utc_first):
            utc_first = data[i].stats.starttime
        if (data[i].stats.endtime > utc_last):
            utc_last = data[i].stats.endtime

    # GPS time of first and last samples (from all blocks)
    gps_first = utc2gps(utc_first)
    gps_last = utc2gps(utc_last)

    if (dt < 1):
        gps_first = round(gps_first/dt)*dt
        gps_last = round(gps_last/dt)*dt
    else:
        gps_first = round(gps_first)
        gps_last = round(gps_last)

    # Check if there is data which falls between the specified GPS times.
    if (gps_last < day_GPSstart or gps_first > day_GPSend):
        # If not, return.
        return chanObj
    else:
        # Set data_flag to True.
        chanObj.data_flag = True

    # We want to figure out the time t_start of the first sample for the day.
    # Example: the first actual data sample (for 1 Hz data) is at gps_first=1.5,
    # This indicates that we should put a zero (for empty data) at t=0.5 and
    # t_start should be 0.5.
    if (gps_first - dt >= day_GPSstart):
        t_start = day_GPSstart + ((gps_first - day_GPSstart) % dt)
    elif (gps_first < day_GPSstart):
        t_start = gps_first + math.ceil((day_GPSstart - gps_first)/float(dt))*dt
    else:
        t_start = gps_first

    print("{0:.8f} {1:.8f}".format(gps_first,t_start))

    # Define time array.
    t_end = t_start + dt*(nSamples - 1) # time of last sample
    t_array = np.linspace(t_start,t_end,num=nSamples,endpoint=True)

    # Loop over blocks and fill the data array.
    for i in range(0,len(data)):
       
        # Start/end times from data block rounded to nearest second or dt,
        # depending on sampling frequency.
        block_sGPS = utc2gps(data[i].stats.starttime)

        if (dt < 1):
            block_sGPS = round(block_sGPS/dt)*dt
        else:
            block_sGPS = round(block_sGPS)
        block_eGPS = block_sGPS + dt*(data[i].stats.npts-1)

        # Code for handling days where a leap second occurs.
        # The Q330/Antelope seems to implement a break (i.e. a new block)
        # near the time of the leap second (either before or after).
        # Before case:
        if (ls_today and (day_GPSend - block_sGPS) < 30 and i > 0):
            if (data[i].stats.starttime == data[i-1].stats.endtime):
                block_sGPS += 1; block_eGPS += 1

        # After case:
        if (ls_yest and (block_eGPS - day_GPSstart) < 30 and i == 0):
            if (data[i].stats.endtime == data[i+1].stats.starttime or
                data[i].stats.endtime == data[i+1].stats.starttime - chanObj.dt + 1):
                block_sGPS -= 1; block_eGPS -= 1

        # The block may contain data which does not fall within the correct GPS times.
        # We check this - d_idx1 and d_idx2 are the first and last indices of data
        # that is in bounds.  data_sGPS and data_eGPS are the times of the first
        # and last sample in the block that fall in bounds.
        if (block_sGPS < day_GPSstart):
            d_idx1 = int((day_GPSstart - block_sGPS)/dt)
            data_sGPS = day_GPSstart
        else:
            d_idx1 = 0
            data_sGPS = block_sGPS

        if (block_eGPS >= day_GPSend):
            d_idx2 = int(data[i].stats.npts - (block_eGPS - (day_GPSend-dt))/dt - 1)
            data_eGPS = day_GPSend-dt
        else:
            d_idx2 = int(data[i].stats.npts - 1)
            data_eGPS = block_eGPS

        # Get indices of the chanObj.data array where we should put the data by comparing
        # data times to the times in t_array, which specifies sample times of chanObj.data.
        # For channels with dt > 1 (i.e., Q or V channels), different blocks in a file can
        # get off relative to each other.  Example: block 1 has samples at 0, 20, 40, etc.
        # and block 2 has samples at 5, 25, 45, etc.
        # This isn't a big deal since the Q** and V**  channels are not very important,
        # so we just adjust the sample times.
        if (dt > 1):
            idx1 = np.where(t_array == data_sGPS)[0]
            if (len(idx1) == 0):
                # Round.
                temp = t_array-data_sGPS
                idx1 = np.argmax(temp[temp < 0])
            else:
                idx1 = int(idx1[0])
        else:
            #print chanObj.fname, i, data_sGPS, (data_sGPS-1107129616.0), t_array[0]-1107129616.0
            idx1 = int(np.where(t_array == data_sGPS)[0][0])
        
        # Get end index for chanObj.data array.
        idx2 = idx1 + (d_idx2 - d_idx1)

        # Copy data from block into channel object data array.
        chanObj.data[idx1:(idx2+1)] = data[i].data[d_idx1:(d_idx2+1)]

    # For days with a leap second, we may need to pull the first sample
    # from the next day's data.
    if (ls_today and data_eGPS < (day_GPSend - dt) and
        (chanObj.channel.find('H',0) == 0 or chanObj.channel.find('L',0) == 0)):
        # Get first second of data from next day and add it to end of chanObj.data.
        addLeapSecondDataFromTomorrow(chanObj,fname,day_UTCstart)

    # Add stuff to chanObj.
    chanObj.startTime = t_start

    # Return.
    return chanObj 

#################################################################################
# Takes data from a miniSEED file (fname) and copies it to a chronologically
# ordered array with uniform sampling time (chanObj.data). Can be pretty tricky
# to deal with multiple blocks in a miniSEED file, leap seconds, etc., so this
# function is a bit complicated.
# INPUT:
#  chanObj        - object holding channel metadata and will have data copied into it.
#  fname          - miniSEED file name.
#  day_UTCstart   - start time of the day for which we are loading data as a
#                   UTCDateTime object.
#  day_index      - index corresponding to the dayvolume we are loading data for.
#                   day_index = 0 means we are loading the miniSEED file corresponding
#                   to the day we want data for, -1 means we are loading the file for
#                   the previous day, etc.  This is because the dayvolumes sometimes
#                   do not end at exactly the end of the day and may contain data
#                   for other days.
# OUTPUT:
#  chanObj - object which hopefully now contains channel data as well as metadata.
def getDataFromMSEED(chanObj,fname,day_UTCstart,day_index):
    # Imports
    import os, math, datetime, sys, pdb
    import numpy as np
    from obspy.core import read
    from timing_utils import utc2gps, getLeapSeconds

    # Set up shorthand variables.
    dt = chanObj.dt

    # Setup: get GPS times of start and end of day and check if a
    # leap second was implemented today or yesterday.
    day_GPSstart, day_length, ls_today, ls_yest = getDayInfo(day_UTCstart)

    # Load data.
    data = read(fname)

    # Get times of first and last data samples which will be stored in chanObj.data.
    # t_first and t_last are RELATIVE to day_GPSstart. All times after this point
    # will be relative to day_GPSstart to deal with rounding errors.
    t_first, t_last, data_flag = getDataTimeRange(data,dt,day_GPSstart,day_length)

    # Check if there is data which falls between the specified GPS times.
    if (data_flag):
        # If there is, set the channel object's data_flag to True.
        chanObj.data_flag = True
    else:
        # If not, return.
        return chanObj
        
    # Loop over blocks and fill the data array.
    for i in range(0,len(data)):
       
        # Get relative times corresponding to first and last data samples in block.
        block_tfirst, block_tlast = getBlockTimes(data[i].stats,dt,day_GPSstart)

        # Code for adjusting times around a leap second. Based on experience with
        # leap second on 6/30/2015, may need to be adjusted in the future.
        if (ls_today or ls_yest):
            block_tfirst, block_tlast = \
                adjustTimesForLeapSecond(data, i, dt, block_tfirst, block_tlast, 
                                         day_GPSstart, day_length, ls_today, ls_yest)

        # Get indices corresponding to the data in the miniSEED block that we want to
        # copy to chanObj.data. Also get corresponding GPS times for debugging.
        c_idx1, c_idx2, data_tfirst, data_tlast = \
            getBlockIndices(block_tfirst, block_tlast, data[i].stats, dt, 
                            day_GPSstart, day_length)
            
    
        # If getBlockIndices finds that the block contains no relevant data,
        # data_sGPS will be -1.  If it's not, then we proceed.
        if (data_tfirst is not None):

            # Get indices corresponding to the elements of chanObj.data which we will
            # will copy the block data into.
            idx1, idx2 = getDataIndices(data_tfirst,t_first,t_last,dt,c_idx2-c_idx1)

            #print data[i].stats
            #print data_tfirst, data_tlast, t_first, t_last, c_idx1, c_idx2, idx1, idx2
            # Copy data from block into channel object data array.
            chanObj.data[idx1:(idx2+1)] = data[i].data[c_idx1:(c_idx2+1)]

    # For days with a leap second, we may need to pull the first sample
    # from the next day's data.  We only do this for H and L data, the other
    # channels are not particularly important.
    if (ls_today and data_tlast < (day_length - dt) and day_index == 0 and
        (chanObj.channel.find('H',0) == 0 or chanObj.channel.find('L',0) == 0)):
        
        # Get data from next day and add it to the end of chanObj.data.
        print "Adding leap second data from tomorrow."
        addLeapSecondDataFromTomorrow(chanObj,fname,day_UTCstart)

    # Add time of first sample to chanObj.
    chanObj.startTime = t_first + day_GPSstart

    # Return.
    return chanObj 


# EOF
