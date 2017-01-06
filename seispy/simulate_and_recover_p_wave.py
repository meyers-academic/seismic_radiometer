from __future__ import division
from matplotlib import use,rc
use('agg')
rc('text',usetex=True)
import numpy as np
from seispy.utils import orf_p_directional
from seispy.simulate import simulate_pwave
from scipy.sparse.linalg import lsqr
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
verbose = False
def set_vec(letter):
    if letter == 'E':
        v = [1,0,0]
    elif letter=='N':
        v = [0,1,0]
    else:
        v = [0,0,1]
    return v

def main():
    conf=0.90
    letters = 'ENZ'
    phi = np.pi / 3
    theta = np.pi / 3
    nstations = 25
    stations = {}
#    stations = {0:[235.6,225.6,225.6],
#	1:[225.7, 297.8, 135.0],
#	2:[537.5, 983.3, 439.6],
#	3:[989.1, 89.2 ,175.5],
#	4:[897.0, 728.6, 950.1],
#	5:[816.3, 891.4, 231.3],
#	6:[151.4, 520.9, 708.4],
#	7:[126.4, 503.7, 812.4]};
    # spiral orientation for stations
    for ii in range(nstations):
        stations[ii] = 1000 * np.array([np.cos(0.2*np.pi*ii),
            np.sin(0.2*np.pi*ii), ii/5])

    for key in stations.keys():
        stations[key]=np.array(stations[key])

    print 'Getting simulated data...'
    data = simulate_pwave(stations, 10, phi, theta, 1, 2000, noise_amp=20, c=5700)
    First = 1
    print 'Populating matrices...'
    counter = 0
    for ii in range(nstations):
        for jj in range(nstations):
            for kk in range(3):
                for ll in range(3):
                    if int(jj)<int(ii):
                        continue
                    elif ll<kk:
                        continue
                    else:
                        p1 = data[ii][letters[kk]].average_fft(fftlength=2)
                        p2 = data[jj][letters[ll]].average_fft(fftlength=2)
                        idx = np.where(p1.frequencies.value==1)
                        p12 = np.conj(p1[idx])*p2[idx]
                        p1 = 1
                        p2 = 1
                    gamma, phis, thetas =\
                        orf_p_directional(set_vec(letters[kk]),set_vec(letters[ll]),stations[jj],
                            stations[ii], 5700, 1)
                    gamma_shape = gamma.shape
                    gamma = gamma.reshape((gamma.size,1))
                    if jj>ii and counter < 1:
#                        print np.asarray(stations[jj])-np.asarray(stations[ii])
#                        print gamma[:5]
                        counter += 1
                    if First == 1:
                        GG = np.dot(np.conj(gamma), np.transpose(gamma))
                        GY = np.conj(gamma)*p12
                        First = 0
                    else:
                        GG += np.dot(np.conj(gamma), np.transpose(gamma))
                        GY += np.conj(gamma)*p12

    print 'GG:'
    print GG[:5,:5]
    print 'GY:'
    print GY[:5]
    print GG.shape
    print GY.shape
    print 'Solving for maximum Likelihood...'
    dphi = phis[2] - phis[1]
    dtheta = thetas[2]-thetas[1]
    S = lsqr(np.real(GG), np.real(GY), iter_lim=1000, atol=1e-6, btol=1e-6)
    print 'Stopped at iteration number ' + str(S[2])
    if S[1]==1:
        print "We've found an exact solution"
    if S[1]==2:
        print "We found an approximate solution"
    print 'Converged to a relative residual of '+str(S[3] /
            np.sqrt((np.abs(GY)**2).sum()))
    S2 = np.copy(S[0])
    S2[S2<0] = 0
    args = np.argsort(S2)[::-1]
    cdf = S2[args].cumsum()
    cdf = cdf/cdf[-1]
    cdf_map = np.zeros(cdf.size)
    cdf_map[args]=cdf
    conf_args = args[cdf<conf]
    S_conf_percent = np.zeros(cdf.size)
    S_conf_percent[conf_args]=1
    S = S[0].reshape(gamma_shape)
    S_conf = S_conf_percent.reshape(gamma_shape)
    print 'Total power in %4.2f confidence region: ' % conf
    print S2[conf_args].sum()
    print S2.size
    print np.sum(cdf<conf)/S2.size

    print 'Plotting...'
    # plot map
    plt.figure()
    plt.subplot(111,projection='aitoff')
    plt.pcolormesh(phis-np.pi-dphi/2, thetas-np.pi/2+dtheta/2, S.T, cmap='viridis')
    plt.colorbar(label='power!')
    CS = plt.contour(phis-np.pi-dphi/2, thetas-np.pi/2+dtheta/2,
            S_conf.T,
            colors='k', linewidth=4,levels=[0])
    plt.scatter(phi, np.pi/2 - theta, s=64, alpha=0.5)
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\theta$')
    plt.grid(True)
    plt.savefig('test_map')
    plt.close()
    # plot phi and theta pdfs
    pdf_map = S2.reshape(gamma_shape) / S2.sum()
    phi_pdf = np.sum(pdf_map, axis=1)
    theta_pdf = np.sum(pdf_map, axis=0)
    plt.figure()
    plt.plot((phis-np.pi)*180/np.pi, phi_pdf)
    plt.plot([phi * 180/np.pi,phi * 180/np.pi], [0,phi_pdf.max()], color='r',
            alpha=0.5, linewidth=4)
    plt.xlabel(r'$\phi$')
    plt.ylabel('PDF')
    plt.savefig('phi_pdf')
    plt.close()
    plt.figure()
    plt.plot((thetas-np.pi/2)*180/np.pi, theta_pdf)
    plt.plot([90- (theta * 180/np.pi),90 - (theta * 180/np.pi)], [0,theta_pdf.max()], color='r',
            alpha=0.5, linewidth=4)
    plt.xlabel(r'$\theta$')
    plt.ylabel('PDF')
    plt.savefig('theta_pdf')
    plt.close()
    # plot only regions within confidence interval
    plt.figure()
    plt.subplot(111,projection='aitoff')
    plt.pcolormesh(phis-np.pi-dphi/2, thetas-np.pi/2+dtheta/2, S_conf.T, cmap='viridis')
    plt.colorbar(label='power!')
    plt.scatter(phi, np.pi/2 - theta, s=64, alpha=0.5)
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\theta$')
    plt.grid(True)
    plt.savefig('test_map_colors')
    plt.close()
    # plot cdf for pixels
    plt.figure()
    plt.plot(cdf)
    plt.plot([0, cdf.size],[0.90, 0.90])
    plt.savefig('cdf')
    # plot station array on 3d plot
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    for key in stations.keys():
        ax.scatter(stations[key][0],stations[key][1],stations[key][2])
    plt.savefig('station_array')
    print 'Done.'



if __name__=='__main__':
    main()
