#! /usr/bin/env python
# Imports
import numpy as np

# Python object for holding data from a Q330 channel and
# the relevant metadata.
# Usage: dataObj = dataChannel(sta_name,chan_name,filename,nSamples)
#        dataObj = dataChannel(sta_name,chan_name,filename,nSamples,startTime,dt,data,data_flag)
class dataChannel:
    
    # Channel count and class description.
    chanCount = 0
    'Class for data from a seismic station and channel.'
    
    def __init__(self,station,channel,nSamples,startTime=-1,dt=None,data=None,data_flag=False):
        self.station = station # station name
        self.channel = channel # chan name
        self.startTime = startTime # GPS time of first sample.
        self.nSamples = nSamples # number of samples.
        self.dt = dt # in seconds
        if (data==None):
            # Fill with zeros
            self.data = np.zeros(nSamples)
        else:
            self.data = data # data (numpy array)
        self.data_flag = data_flag # 0 = no actual data, array contains filler values, 1 = actual data present
        dataChannel.chanCount += 1

    def displayInfo(self):
        GPS_end = self.startTime + self.dt*(self.nSamples-1)
        if (self.data_flag == 1):
            print "Station " + self.station + ", channel " + \
                self.channel + ".  Data starts at GPS time " + \
                str(self.startTime) + " and has " + str(self.data.size) + \
                " samples with dt " + str(self.dt) + " seconds."
        else:
            print "No data found for station " + self.station + ", channel " + \
                self.channel + " between GPS times " + str(self.startTime) + " and " \
                + str(GPS_end) + " (inclusive).  The data array contains filler values, and has " + \
                str(self.data.size) + " samples with dt " + str(self.dt) + " seconds."

    def displayCount(self):
        print "Total number of dataChannels: " + str(dataChannel.chanCount)

# EOF
