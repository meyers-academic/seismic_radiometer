import ConfigParser
from collections import OrderedDict
#from termcolor import colored

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def config_section_map(c, section):
    dict1 = OrderedDict()
    options = c.options(section)
    for option in options:
        dict1[option] = c.get(section, option).strip()
    return dict1


def config_list(c):
    dict1 = {}
    for section in c.sections():
        dict1[section] = config_section_map(c, section)
    return dict1


def config_pipeline_ini(c):
    dict1 = OrderedDict()
    for section in c.sections():
        dict1[section] = config_section_map(c, section)
    return dict1


def read_ini(file):
    """
    Read ini file for analysis and put
    it into a dict.

    Parameters
    ----------
    file : str
        ini file for coherence pipeline

    Returns
    -------
    dict1 : dict
        dictionary with parameters for pipeline
    """
    c = ConfigParser.ConfigParser()
    c.optionxform=str
    c.read(file)
    dict1 = config_pipeline_ini(c)
    return dict1

def print_params(params):
    for key in params.keys():
        print '\n%s PARAMS' % key.upper()
        for ii,param in enumerate(params[key].keys()):
            if ii % 2 ==0:
                print('\t{:<20}{:<20}'.format(*[param,
                        params[key][param]]))
            else:
                print('\t{:<20}{:<20}'.format(*[param,
                        params[key][param]]))

def try_and_set(params, section, parameter, default_val=None, warning=False,
        special=None):
    try:
        params[section][parameter]
        if special=='comma_list':
            val = params[section][parameter]
            params[section][parameter] = [float(v) for v in val.split(',')]
    except KeyError:
        if default_val is None:
            raise KeyError('Must set %s in config file' % parameter)
        else:
            params[section][parameter] = default_val
    if warning:
        print 'WARNING: setting %s to %s' % (parameter, str(default_val))
    return params
