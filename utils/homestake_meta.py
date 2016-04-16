#! /usr/bin/env python

# List of stations.
# type can be 'all': all stations, including MAST.
#             'surf': surface stations only.
#             'ug': underground stations only.
def getStations(type):

    # stations
    sta_ug = ['300','800','1700','A2000','B2000','C2000','D2000','E2000','A4100',
              'C4100','D4100','A4850','B4850','C4850','D4850']
    sta_surf = ['DEAD','LHS','ORO','ROSS','RRDG','SHL','TPK','YATES','WTP']
    sta_other = ['MAST']

    if (type.lower() == 'all'):
        sta_list = sta_ug + sta_surf + sta_other
    elif ((type.lower() == 'ug') or (type.lower() == 'underground')):
        sta_list = sta_ug
    elif ((type.lower() == 'surf') or (type.lower() == 'surface')):
        sta_list = sta_surf
    else:
        raise ValueError('type should be \'all\' or \'ug\' or \'surf\'.')

    return sta_list


# List of all channels from Q330 and Antelope.
def getChannels(type):

    # List of all channels saved by Antelope
    # HHZ - Channel 1 (z-axis) (100 Hz)
    # HHE - Channel 2 (East) (100 Hz)
    # HHN - Channel 3 (North) (100 Hz)
    # LHZ - Channel 1 (z-axis) (1 Hz)
    # LHE - Channel 2 (East) (1 Hz)
    # LHN - Channel 3 (North) (1 Hz)
    # LCQ - internal clock quality percentage (0-100)
    # LCE - clock phase error (microseconds)
    # LCC - GPS clock quality
    # LPL - clock phase lock loop status
    # LCL - time since GPS lock was lost
    # QBD - total number of Q330 reboots in last 24 hours
    # QBP - logical port buffer percent full from real-time status
    # QDG - data gaps (seconds)
    # QDL - current data latency (seconds)
    # QDR - current total input+output data rate (bits/second)
    # QEF - overall communications efficiency (percent)
    # QG1 - total number of data gaps in last 1 hour
    # QGD - total number of data gaps in last 24 hours
    # QID - total number of datalogger ip-address cahnges in last 24 hours
    # QLD - total number of comm link cycles in last 24 hours
    # QPD - total number of POCs received in last 24 hours
    # QRD - total number of bytes read in last 24 hours
    # QRT - current run time (seconds)
    # QTH - current throttle setting (bits/second0
    # QTP - ratio of seconds read to real-time clock
    # QWD - total number of bytes written in last 24 hours
    # VCO - voltage controlled oscillator value
    # VEA - antenna current
    # VEC - main system current
    # VEP - main system voltage
    # VKI - main system temperature
    # VM1 - mass position for channel 1
    # VM2 - mass position for channel 2
    # VM3 - mass position for channel 3
    # VPB - percentage of packet buffer full
    # VTW - main system opto inputs
    # note: all V** channels are updated every 10 seconds

    data_chan = ['HHE','HHN','HHZ','LHE','LHN','LHZ']
    status_chan = ['LCE','LCQ','LPL','LCL','LCC','QBD','QBP','QDG','QDL',
                   'QDR','QEF','QG1','QGD','QID','QLD','QPD','QRD','QRT',
                   'QTH','QTP','QWD','VCO','VEA','VEC','VEP','VKI','VM1',
                   'VM2','VM3','VPB','VTW']
    useful_chan = ['HHE','HHN','HHZ','LCE','LCQ','VEP','VKI','VM1','VM2','VM3']

    if (type.lower() == 'all'):
        chan_list = data_chan + status_chan
    elif (type.lower() == 'data'):
        chan_list = data_chan
    elif (type.lower() == 'status'):
        chan_list = status_chan
    elif (type.lower() == 'useful'):
        chan_list = useful_chan
    else:
        raise ValueError('type should be \'all\' or \'data\' or \'status\' or \'useful\'.')

    return chan_list

# EOF
