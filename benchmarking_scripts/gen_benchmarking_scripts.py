

types = ['p', 's', 'r', 'ps', 'pr', 'sr', 'psr']
velocities = ['5700', '2500', '3000',
              '5700,3000','5700,2500',
              '3000,2500','5700,3000,2500']
phis = {}
phis = {'p':120, 's':240, 'r':180}
thetas = {'p':120, 's':60, 'r':90}
vels = {'p':5700, 's':3000, 'r':2500}

recovery_section = """
[Recovery]
recovery_string={rec_str}
velocities={vel_str}
nproc=8
array_type=homestake
tag={tag}
frequency=1
output_directory=./benchmarking_scripts/{dir}
"""

injection_section = """
[Injection {inj_num}]
type={inj_type}
phi={phi}
theta={theta}
frequency=1
velocity={vel}
amplitude=1e-4
"""

f2 = open('run_benchmarking.sh', 'w')
f2.write('#!/bin/bash\n')
for ii,rec_str in enumerate(types):
    for rec_str2 in types:
        fname = '%s_rec_%s_inj.ini' % (rec_str, rec_str2)
        f = open(fname, 'w')
        recovery_kwargs = {}
        recovery_kwargs['rec_str'] = rec_str
        recovery_kwargs['vel_str'] = velocities[ii]
        recovery_kwargs['tag'] = '%s_rec_%s_inj' % (rec_str, rec_str2)
        recovery_kwargs['dir'] = '%s_rec_%s_inj' % (rec_str, rec_str2)
        f.write(recovery_section.format(**recovery_kwargs))
        for jj,string in enumerate(rec_str2):
            inj_kwargs = {}
            inj_kwargs['inj_type'] = string
            inj_kwargs['inj_num'] = str(jj+1)
            inj_kwargs['phi'] = phis[string]
            inj_kwargs['theta'] = thetas[string]
            inj_kwargs['vel'] = vels[string]
            f.write(injection_section.format(**inj_kwargs))
        f.close()
        f2.write("echo '%s recovery %s injected'\n" % (rec_str, rec_str2))
        f2.write("mkdir benchmarking_scripts/%s_rec_%s_inj\n" % (rec_str,
            rec_str2))
        f2.write("python -m seispy.engine -c './benchmarking_scripts/%s' >\
'benchmarking_scripts/%s_rec_%s_inj/%s_rec_%s_inj_summary.txt'\n"
                % (fname, rec_str, rec_str2, rec_str, rec_str2))

f2.close()
