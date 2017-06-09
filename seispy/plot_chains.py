import matplotlib
matplotlib.use('agg')
from corner import corner
import numpy as np

dat = np.loadtxt('chains/1-post_equal_weights.dat')
figure = corner(dat[:, :11], labels=['$\phi_p$', '$\\theta_p$', 'v_p', 'A_p',
    '$\phi_{sh}$', '$\\theta_{sh}$', 'v_s', 'A_{sh}',
    '$\phi_{sv}$','$\\theta_{sv}$', 'A_{sv}'], quantiles =
        [0.16, 0.5, 0.84], show_titles=True)

figure.savefig('mycorner_p_and_s')
