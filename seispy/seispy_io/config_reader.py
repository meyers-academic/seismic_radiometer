from __future__ import division
import ConfigParser
from .utils import *

def read_config(cfile):
    params = read_ini(cfile)
    params = try_and_set(params, 'Recovery','phimesh',10)
    params = try_and_set(params, 'Recovery','thetamesh',10)
    params = try_and_set(params, 'Recovery','segdur',100)
    params = try_and_set(params, 'Recovery','duration',100)
    params = try_and_set(params, 'Recovery','velocities',special='comma_list')
    params = try_and_set(params, 'Recovery','frequency')
    params = try_and_set(params, 'Recovery','recovery_string')
    params = try_and_set(params, 'Recovery','array_type')
    params = try_and_set(params, 'Recovery','output_directory','./')
    params = try_and_set(params, 'Recovery','num_stations', 12)
    params = try_and_set(params, 'Recovery','alpha', 1000)
    params = try_and_set(params, 'Recovery','epsilon', 0.1)
    params = try_and_set(params, 'Recovery','nproc', 1)
    params = try_and_set(params, 'Recovery','tag', 'test')
    params = try_and_set(params, 'Recovery','sample_rate', 100)
    for key in params.keys():
        sp = key.split(' ')
        if sp[0]=='Injection':
            params = try_and_set(params, key, 'phi')
            params = try_and_set(params, key, 'theta')
            params = try_and_set(params, key, 'velocity')
            params = try_and_set(params, key, 'frequency')
            params = try_and_set(params, key, 'amplitude')
            params = try_and_set(params, key, 'alpha', default_val=1000)
            params = try_and_set(params, key, 'epsilon', default_val=0.1)
            params = try_and_set(params, key, 'psi', default_val=0)
            params = try_and_set(params, key, 'phase', default_val=0)
            params = try_and_set(params, key, 'sample_rate', default_val=100)
            params = try_and_set(params, key, 'type')
    return params

