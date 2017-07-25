import numpy as np
import glob
from gwpy.timeseries import TimeSeries

#  set station locations
locs = {}
locs['DEAD'] = [44.382699, -103.7532, 1498, 0]
locs['LHS'] = [44.347299, -103.7748, 1684, 0]
locs['ORO'] = [44.343499, -103.7523, 1543, 0]
locs['ROSS'] = [44.345099, -103.7576, 1628, 0]
locs['RRDG'] = [44.359499, -103.7654, 1677, 0]
locs['SHL'] = [44.316899, -103.7098, 1772, 0]
locs['TPK'] = [44.340799, -103.7980, 1740, 0]
locs['WTP'] = [44.353899, -103.7423, 1555, 0]
locs['YATES'] = [44.352199, -103.7515, 1625, 0]
locs['300'] = [44.346399, -103.7569, 1505.1, 91.44]
locs['800'] = [44.347399, -103.7589, 1350, 243.84]
locs['1700'] = [44.351699, -103.7584, 1073.8, 518.16]
locs['A2000'] = [44.351099, -103.7623, 983.3, 609.6]
locs['B2000'] = [44.349199, -103.7605, 983.3, 609.6]
locs['C2000'] = [44.351599, -103.7637, 983.1, 609.6]
locs['D2000'] = [44.353299, -103.7689, 983.0, 609.6]
locs['E2000'] = [44.356399, -103.7716, 983.1, 609.6]
locs['A4100'] = [44.344899, -103.7535, 342.5, 1249.68]
locs['C4100'] = [44.351199, -103.7507, 342.3, 1249.68]
locs['D4100'] = [44.343599, -103.7511, 342.5, 1249.68]
locs['A4850'] = [44.340499, -103.7625, 115.2, 1478.28]
locs['B4850'] = [44.346299, -103.7581, 114.9, 1478.28]
locs['C4850'] = [44.346299, -103.7528, 114.6, 1478.28]
locs['D4850'] = [44.352999, -103.7505, 115.2, 1478.28]

# read frame
def read_frame(frame, channel, st=None, et=None, cfac=1.589459e-9):
    """
    reads ligo frames

    Parameters
    ----------
    frame : `str`
        filepath to a frame
    channel : `str`
        channel in the frame to load
    st : `int`, date string, optional
        optional start time. defaults to beginning
        of frame
    et : `int ,date string, optional
        optional end time. defaults to end
        of frame

    Returns
    -------
    TS : `Trace`
        time series trace
    """

    if st is not None and et is not None:
        d1 = cfac * TimeSeries.read(frame, channel, st, et).detrend()
    else:
        d1 = cfac * TimeSeries.read(frame, channel).detrend()
    d1.location = d1.get_location()
    return d1

def fetch_homestake_data(st, et, channel, framedir='./'):
    """
    fetch data based on location of frames

    Parameters
    ----------
    st : `int`
        start time (GPS time)
    et : `int`
        end time (GPS time)
    channel : `channel`
        channel to load data for

    Returns
    -------
    TS : `gwpy.timeseries.TimeSeries`
        TimeSeries object containing data between start and end times
    """
    # uncomment when not testing
    # for looping over directories where frames
    # are located
    st_dir = int(str(st)[:5])
    et_dir = int(str(et)[:5])
    dirs = np.arange(st_dir, et_dir + 1)
    files = []
    for directory in dirs:
        # print 'FRAME READING TESTING MODE!!!'
        loaddir = '%s/M-%d/' % (framedir, directory)
        new_files = sorted(glob.glob(loaddir + '/*.gwf'))
        files.extend(new_files)
    vals = np.asarray([])
    if len(files)==0:
        raise ValueError('No files found...we looked here: %s' % loaddir)

#    for file in files:
    for ii in range(len(files)):
        file = files[ii]
        # start is before frame start, end is before frame end
        # we want to load fst -> et
        fst = int(file.split('-')[-2])
        dur = int(file.split('-')[-1][:-4])
        if st <= fst and et <= fst + dur and et >= fst:
            val = read_frame(file, channel, st=fst, et=et)
            vals = np.hstack((vals, val.value))
        # start is after frame start, end is before frame end
        # we want to load only st -> et
        elif st >= fst and et <= fst + dur:
            val = read_frame(file, channel, st=st, et=et)
            vals = np.hstack((vals, val.value))
        # start is after frame start end is after or equal to frame end
        # we want to load st -> fst + dur
        elif st >= fst and st < (fst + dur) and et >= fst + dur:
            val = read_frame(file, channel, st=st, et=fst + dur)
            vals = np.hstack((vals, val.value))
        # start is before frame start, end is after frame end
        # load fst -> fst + dur (whole frame)
        elif st <= fst and et >= fst + dur:
            val = read_frame(file, channel)
            vals = np.hstack((vals, val.value))
        else:
            continue
    TS = TimeSeries(vals, x0=st, dx=val.dx, name=val.name, channel=val.channel)
    loc = TS.get_location()
    TS.location = loc
    return TS

def list_station_names():
    print locs.keys()


def get_coordinates(station):
    """
    Gets homestake station coordinates

    Returns
    -------
    coordinate : `list`
        list [latitude, longitude, elevation, depth]
    """
    staname = station.split(':')[0]
    return np.asarray(locs[staname])


