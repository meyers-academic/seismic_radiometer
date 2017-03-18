"""
Class for seismic event metadata. Includes latitude, longitude, and time of
event, along with an ID number. A user can also define a time window
(relative to the event time) and taper lengths for data processing.

(Originally written by Tanner Prestegard)
"""

# Imports
import datetime
# Check if obspy is available.
try:
    # Import stuff from obspy.
    from obspy.core.utcdatetime import UTCDateTime
except ImportError:
    raise ImportError('Error: can\'t find obspy.  Please install it or add it to your $PYTHONPATH.')

# chanObj class.
class Event:
    """Class for seismic events from an IRIS database or our own database."""

    def __init__(self, latitude, longitude, time, evID=None, magnitude=None,
                 win_start=0, win_end=-1, taper_start=10, taper_end=10):
        # Defining and processing class members.
        self.latitude = float(latitude) # Latitude of event. (numeric)
        self.longitude = float(longitude) # Longitude of event. (numeric)
        self.time = time # time of event (will convert to UTCDateTime if not already)
        self.evID = evID # event ID number.
        self.magnitude = magnitude # event magnitude (for earthquakes)
        self.win_start = win_start # start of data we want to look at (relative to event time).
        self.win_end = win_end # end of data we want to look at (relative to event time).
        self.taper_start = taper_start # length of starting taper (s)
        self.taper_end = taper_end # length of ending taper (s)

        # Make time a UTCDateTime object.
        if isinstance(time, str):
            self.time = UTCDateTime(datetime.datetime.strptime(time, '%m/%d/%Y %H:%M:%S'))
        elif isinstance(time, UTCDateTime):
            self.time = time
        else:
            raise TypeError('time should be a string formatted like MM/DD/YYYY ' +
                            'HH:MM:SS or a UTCDateTime object.')


    def displayInfo(self):
        """Displays info about event, useful for debugging."""
        print(self.__module__ + ":")
        print("\tEvent ID: " + (str(self.evID) if self.evID is not None else 'N/A'))
        print("\tLatitude: " + str(self.latitude))
        print("\tLongitude: " + str(self.longitude))
        print("\tEvent time: " + self.time.strftime('%H:%M:%S.%f %m/%d/%Y'))
        print("\tWindow: " + str(self.win_start) + " - " + str(self.win_end) + " sec.")
        print("\tTaper: " + str(self.taper_start) + " sec. start, " + str(self.taper_end) + " sec. end.")
        if (self.magnitude is None):
            print("\tMagnitude: N/A")
        else:
            print("\tMagnitude: " + str(self.magnitude))

    def createLog(self, filename):
        """Generate log file with event information, used for web interface."""
        # Set up lines for writing to file.
        date_line = "Date " + self.time.strftime('%m/%d/%Y') + "\n"
        time_line = "Time " + self.time.strftime('%H:%M:%S.%f') + "\n"
        lat_line = "Latitude " + str(self.latitude) + "\n"
        long_line = "Longitude " + str(self.longitude) + "\n"
        mag_line = "Magnitude " + (str(self.magnitude) if self.magnitude is not None else "N/A") + "\n"

        # write to file.
        f = open(filename,'w')
        f.write(date_line)
        f.write(time_line)
        f.write(lat_line)
        f.write(long_line)
        f.write(mag_line)
        f.close()

        # End of function


# EOF