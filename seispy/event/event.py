"""
Class for seismic event metadata. Includes latitude, longitude, and time of
event, along with an ID number. A user can also define a time window
(relative to the event time) and taper lengths for data processing.

Includes a built in table from gwpy.table.Table

(Originally written by Tanner Prestegard)
"""

# Imports
import datetime
import os, csv
from gwpy.table import Table
from collections import OrderedDict
import numpy as np
# Check if obspy is available.
try:
    # Import stuff from obspy.
    from obspy.core.utcdatetime import UTCDateTime
except ImportError:
    raise ImportError('Error: can\'t find obspy.  Please install it or add it to your $PYTHONPATH.')
from seispy.station import Seismometer
from scipy.optimize import curve_fit
from scipy.signal import argrelmax

analyzed_names = ['latitude', 'longitude', 'time', 'evID',
                  'magnitude', 'win_start', 'win_end', 'taper_start',
                  'taper_end', 'filter frequency', 'peak amplitude', 
                  'peak time', 'peak time minimum',
                  'peak time maximum', 'velocity', 'bearing', 'distance',
                  'channel', 'analyzed', 'station', 'depth']


# TODO we need some unittests for these classes
# chanObj class.
class Event(dict):
    """Class for seismic events from an IRIS database or our own database.
    """

    def __init__(self, latitude, longitude, time, evID=None,  win_start=0, taper_start=10,
                 win_end=-1, taper_end=10, magnitude=None, analyzed=False, **kwargs):
        # Defining and processing class members.
        super(Event, self).__init__(**kwargs)
        self['latitude'] = float(latitude) # Latitude of event. (numeric)
        self['longitude'] = float(longitude) # Longitude of event. (numeric)
        self['time'] = time # time of event (will convert to UTCDateTime if not already)
        self['evID'] = evID # event ID number.
        self['magnitude'] = magnitude # event magnitude (for earthquakes)
        self['win_start'] = win_start # start of data we want to look at (relative to event time).
        self['win_end'] = win_end # end of data we want to look at (relative to event time).
        self['taper_start'] = taper_start # length of starting taper (s)
        self['taper_end'] = taper_end # length of ending taper (s)
        self['analyzed'] = analyzed

        # Make time a UTCDateTime object.
        if isinstance(time, str):
            self['time'] = UTCDateTime(datetime.datetime.strptime(time, '%m/%d/%Y %H:%M:%S'))
        elif isinstance(time, UTCDateTime):
            self['time'] = time
        else:
            raise TypeError('time should be a string formatted like MM/DD/YYYY ' +
                            'HH:MM:SS or a UTCDateTime object.')


    def __repr__(self):
        string =\
        """
        Event ID: {evID}
        Latitude: {latitude}
        Longitude: {longitude}
        Event Time: {time}
        Window: {win_start}-{win_end}
        Taper: {taper_start} sec. start, {taper_end} sec. end.
        Magnitude: {magnitude}
        """.format(**self)
        return string

    @property
    def latitude(self):
        return self['latitude']

    @property
    def longitude(self):
        return self['longitude']

    @property
    def analyzed(self):
        return self['analyzed']

    @property
    def time(self):
        return self['time']

    @property
    def win_start(self):
        return self['win_start']

    @property
    def win_end(self):
        return self['win_end']

    @property
    def taper_start(self):
        return self['taper_start']

    @property
    def taper_end(self):
        return self['taper_end']

    @property
    def magnitude(self):
        return self['magnitude']

    def createLog(self, filename):
        """Generate log file with event information, used for web interface."""
        # Set up lines for writing to file.
        date_line = "Date " + self['time'].strftime('%m/%d/%Y') + "\n"
        time_line = "Time " + self['time'].strftime('%H:%M:%S.%f') + "\n"
        lat_line = "Latitude " + str(self['latitude']) + "\n"
        long_line = "Longitude " + str(self['longitude']) + "\n"
        mag_line = "Magnitude " + (str(self['magnitude']) if self['magnitude'] is not None else "N/A") + "\n"

        # write to file.
        f = open(filename,'w')
        f.write(date_line)
        f.write(time_line)
        f.write(lat_line)
        f.write(long_line)
        f.write(mag_line)
        f.close()

    def analyze(self, station, frequencies, framedir='./',
            return_envelopes=False, changed_bearing=None):
        if isinstance(self.time, UTCDateTime):
            self['time'] = utc2gps(self.time)
        data = Seismometer.fetch_data(station, self.time+self.win_start, self.time+self.win_end,
                                      framedir=framedir, chans_type='fast_chans')
        analyzed_event = self.copy()
        analyzed_event['station'] = station
        # detrend the data
        for key in data.keys():
            data[key] = data[key].detrend()
        # get bearing and distance
        dist, bearing = data['HHZ'].get_wave_path(self)
        analyzed_event['distance'] = dist
        analyzed_event['bearing'] = bearing
        # add transverse and radial channels to this seismometer
        if changed_bearing is not None:
            print 'CHANGING BEARING: %4.2f to %4.2f' % (bearing,
                    changed_bearing)
            data.rotate_RT(changed_bearing)
        else:
            data.rotate_RT(bearing)
        # taper around the rayleigh-wave part
        for key in data.keys():
            data[key] = data[key].taper(event=self)
        final_table = Table(names=analyzed_names)
        final_table['channel'] = final_table['channel'].astype('str')
        final_table['station'] = final_table['channel'].astype('str')
        env_dict = OrderedDict()
        filtered_final = OrderedDict()
        for frequency in frequencies:
            filtered_final[frequency] = {}
            env_dict[frequency] = {}
            for key in data.keys():
                # channel
                analyzed_event['channel'] = key
                # filter data
                filtered = data[key].gaussian_filter(frequency, 0.1)
                # get depth
                analyzed_event['depth'] = data[key].get_coordinates()[-1]
                # add filtered data to dict
                filtered_final[frequency][key] = filtered
                # get hilbert
                env = filtered.hilbert()
                # add envelope to dict
                env_dict[frequency][key] = env
                # get local max indices
                local_max_idxs = argrelmax(env)
                # get global max index
                global_max_idx = np.argmax(env)
                # normalize local maxima by global maximum
                if sum(np.zeros(filtered.size) == filtered.value):
                    raise ValueError('Station %s has all zeros during this\
                    time.' % station)
                try:
                    lm_normed = env[local_max_idxs[0]] / env[global_max_idx]
                except:
                    print local_max_idxs, global_max_idx
                    print env
                    print filtered
                    raise('ValueError')
                # get index of first local maxima that is 90% of global maximum
                try:
                    idx_arrival = np.where(lm_normed > 0.9)[0][0]
                except:
                    idx_arrival = np.argmax(lm_normed)
                # set peak amplitude to that
                analyzed_event['peak amplitude'] =\
                    env[local_max_idxs[0][idx_arrival]]
                # set peak time
                pt = env.times.value[idx_arrival]
                analyzed_event['peak time'] = pt
                # get group velocity based on event time
                analyzed_event['velocity'] = (dist / (pt - self.time))
                # set some bookkeeping variables
                analyzed_event['analyzed'] = True
                analyzed_event['filter frequency'] = frequency
                # add this event to our table
                final_table.add_row(analyzed_event)
        if return_envelopes:
            return final_table, env_dict, filtered_final
        else:
            return final_table


def getEstimateAndConfLevel(times, envelope, conf):
    """
    Get our estimated value and confidence interval given a distribution.
    In this case we're estimating the peak time in this region.

    Parameters
    ----------
    times : `numpy.ndarray`
        possible peak times
    envelope : `numpy.ndarray`
        ampltiude of wave train
    conf : `float`
        confidence level
    Returns
    -------
    conf : `numpy.ndarray`
        All of the values from `times` in our confidence
        interval (or less than our upper limit)
    estimate : `float`
       Our measured estimate of the value of :math:`h_0` from the posterior.
       Nominally hardcoded to be the median of the returned interval.
    """
    # get indices to sort descending
    sort_idxs = np.argsort(envelope)[::-1]
    envelope_sort = envelope[sort_idxs]
    times_sort = times[sort_idxs]
    conf = times_sort[np.where(np.cumsum(envelope_sort) < conf)]
    return conf, np.median(conf)

class EventTable(Table):
    """
    Create an event table
    """
    def __init__(self, **kwargs):
        # Let's add column names to begin with
        super(EventTable, self).__init__(**kwargs)

    @classmethod
    def read_iris_db(cls, db_file, window_file=None):
        ev = Table.read(db_file, format='ascii')
        ev.rename_column('EventID','evID')
        if window_file is not None:
            evIDs = []
            if (os.path.isfile(window_file)):
                wf_raw = load_file(window_file,"|","#")

                # Turn wf into a dictionary.
                wf_tab = Table(names=['evID','win_start',
                'taper_start','win_end','taper_end'])
                for wf_i in wf_raw:
                    if (wf_i[0] not in wf_tab['evID']):
                        evIDs.append(int(wf_i[0]))
                        wf_tab.add_row(wf_i)
                    else:
                        raise KeyError('Duplicate event: key already in wf.')
            else:
                wf_tab = None
        else:
            wf_tab = None
        short_table = Table()
        if wf_tab is not None:
            good_events = evIDs
            long_table = ev[ev['evID']==good_events]
            short_table.add_columns([long_table['Latitude'],long_table['Longitude'],long_table['Time']])
            short_table.add_columns(wf_tab.columns.values())
        else:
            short_table.add_columns([ev['Latitude'],ev['Longitude'],ev['Time']])
        short_table['Time'] = [UTCDateTime(short_table['Time'][ii]) for ii in
                range(len(short_table['Time']))]
        short_table.rename_column('Time','time')
        short_table.rename_column('Latitude','latitude')
        short_table.rename_column('Longitude','longitude')
        return short_table

    @classmethod
    def read_seis_db(cls, db_file, window_file=None):
        ev = Table.read(db_file, format='ascii')
        ev.rename_column('col1','latitude')
        ev.rename_column('col2','longitude')
        ev.rename_column('col4','time')
        # clean up garbage events
        remove_rows = []
        for ii in range(len(ev)):
            if ev['col23'][ii]=='-':
                remove_rows.append(ii)
        print remove_rows
        ev.remove_rows(remove_rows)
        # sort by time
        ev.sort('time')
        # add window information
        if window_file is not None:
            evIDs = []
            if (os.path.isfile(window_file)):
                wf_raw = load_file(window_file,"|","#")

                # Turn wf into a dictionary.
                wf_tab = Table(names=['evID','win_start',
                'taper_start','win_end','taper_end'])
                for wf_i in wf_raw:
                    if (wf_i[0] not in wf_tab['evID']):
                        evIDs.append(int(wf_i[0]))
                        wf_tab.add_row(wf_i)
                    else:
                        raise KeyError('Duplicate event: key already in wf.')
            else:
                wf_tab = None
        else:
            wf_tab = None
        short_table = Table()
        if wf_tab is not None:
            good_events = evIDs
            long_table = ev[good_events]
            short_table.add_columns([long_table['latitude'],long_table['longitude'],long_table['time']])
            short_table.add_columns(wf_tab.columns.values())
        else:
            short_table.add_columns([ev['latitude'],ev['longitude'],ev['time']])
        short_table['time'] = [UTCDateTime(short_table['time'][ii]) for ii in
                range(len(short_table['time']))]
        return short_table, long_table
def load_file(fname, delim, comment_str):
    '''
    Helper function for loading delimited text files.
    Input:
      fname       - filename
      delim       - delimiter (string)
      comment_str - string which denotes a comment (lines begining with
                    this string will be skipped).

    Output:
      db - list of lists, each sublist corresponds to a row.
    '''

    db = []
    with open(fname) as f:
        reader = csv.reader(f, delimiter=delim)
        for row in reader:
            # Read all non-comment lines.
            if (row[0][0] != comment_str):
                db.append(row)

    return db
# Miscellaneous functions for writing frames from miniseed files.

# Takes UTC time (input as an obspy.core.utcdatetime.UTCDateTime object)
# and determines the number of leap seconds that have occurred since
# GPS time began.
# Up-to-date as of July 2015.  Once a new leap second is added, this
# code will need to be updated.  The leap second can be added to the
# list before it is actually implemented and the code should handle
# it without any trouble.
def getLeapSeconds(UTC_time):

    from obspy.core.utcdatetime import UTCDateTime

    # List of leap seconds since GPS zero time (00:00:00, Jan. 6, 1980).
    # References: tf.nist.gov/pubs/bulletin/leapsecond.htm
    #             en.wikipedia.org/wiki/Leap_second
    ls_array = [UTCDateTime(1981,06,30,23,59,59),
                UTCDateTime(1982,06,30,23,59,59),
                UTCDateTime(1983,06,30,23,59,59),
                UTCDateTime(1985,06,30,23,59,59),
                UTCDateTime(1987,12,31,23,59,59),
                UTCDateTime(1989,12,31,23,59,59),
                UTCDateTime(1990,12,31,23,59,59),
                UTCDateTime(1992,06,30,23,59,59),
                UTCDateTime(1993,06,30,23,59,59),
                UTCDateTime(1994,06,30,23,59,59),
                UTCDateTime(1995,12,31,23,59,59),
                UTCDateTime(1997,06,30,23,59,59),
                UTCDateTime(1998,12,31,23,59,59),
                UTCDateTime(2005,12,31,23,59,59),
                UTCDateTime(2008,12,31,23,59,59),
                UTCDateTime(2012,06,30,23,59,59),
                UTCDateTime(2015,06,30,23,59,59),
                UTCDateTime(2016,12,31,23,59,59)]

    # Count how many leap seconds have occurred before UTC_time.
    ls  = [date for date in ls_array if UTC_time > date]

    # Check if a leap second is being added on the data specified
    # by UTC_time.
    ls_today = [date for date in ls_array if \
                    UTC_time.year == date.year and \
                    UTC_time.month == date.month and \
                    UTC_time.day == date.day]

    # Return number of leap seconds.
    return (len(ls), len(ls_today))


# Converts from UTC time (input as an obspy.core.utcdatetime.UTCDateTime
# object) to GPS time (in seconds).
def utc2gps(UTC_time):
    from obspy.core.utcdatetime import UTCDateTime

    # UTC time when GPS time = 0 (start of GPS time)
    time_zero = UTCDateTime(1980,1,6,0)

    # Get total number of UTC seconds since the beginning of GPS time.
    if isinstance(UTC_time, UTCDateTime):
        UTC_seconds = UTC_time.timestamp - time_zero.timestamp

    # Raise error if time requested is before beginning of GPS time.
    if (UTC_seconds < 0):
        err_msg = ('You have specified ' + UTC_time.ctime() +
                   ', which occurred before GPS begin (Jan. 6, 1980).')
        raise ValueError(err_msg);

    # Account for leap seconds to get total GPS seconds.
    num_ls = getLeapSeconds(UTC_time)[0]
    GPS_seconds = UTC_seconds + num_ls

    return GPS_seconds

# End of file
