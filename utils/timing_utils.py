#! /usr/bin/env python

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
    import datetime
    
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
                UTCDateTime(2015,06,30,23,59,59)]

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
