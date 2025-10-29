import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import zHead as hd
from astropy.cosmology import FlatLambdaCDM

arg_options = hd.handle_arguments()
ActiveRunIDs = arg_options[0]
print(ActiveRunIDs)
date = arg_options[1]

###

CALCULATE = False

MODE = hd.DEFAULT_MODE

FSNAP = 100
NSNAPS = 401

gammas = np.loadtxt('data/gamma_dat')[:,1][:6]
gamma_errs = np.loadtxt('data/gamma_dat')[:,2][:6]

def disp( v, m ):
    Mtot = sum(m)
    v_avg = sum( m*v ) / Mtot
    v2_avg = sum( m*v*v ) / Mtot
    disp2 = v2_avg - v_avg*v_avg
    return np.sqrt(disp2)

Nruns = len(hd.runlabels)

if CALCULATE:

    vz_disps_allruns = []
    vr_disps_allruns = []
    for j,(RunName, SnapDir, FofDir) in enumerate(zip(hd.runlabels,hd.snap_dirs,hd.fof_dirs)):
        if j in ActiveRunIDs:

            times = []
            vz_disps = []
            vr_disps = []
            for i in range(FSNAP,NSNAPS):

                print(i)

                f = h5py.File(SnapDir + hd.get_istr3(i) + '.hdf5')

                SimTime = f[u'Header'].attrs[u'Time'] * 1000.

                xyz4 = np.array(f[u'PartType4'][u'Coordinates'])
                m4 = 1.e10 * np.array(f[u'PartType4'][u'Masses'])

                com_x = np.sum(m4*xyz4[:,0]) / np.sum(m4)
                com_y = np.sum(m4*xyz4[:,1]) / np.sum(m4)
                com_z = np.sum(m4*xyz4[:,2]) / np.sum(m4)
                com = np.array([ com_x, com_y, com_z ])
                xyz4 -= com

                r4 = np.sqrt( sum([ xyz4[:,b]*xyz4[:,b] for b in range(3) ]) )
                r4_Cyl = np.sqrt( sum([ xyz4[:,b]*xyz4[:,b] for b in range(2) ]) )

                vel4 = np.array(f[u'PartType4'][u'Velocities'])

                com_vx = np.sum(m4*vel4[:,0]) / np.sum(m4)
                com_vy = np.sum(m4*vel4[:,1]) / np.sum(m4)
                com_vz = np.sum(m4*vel4[:,2]) / np.sum(m4)
                v_cen = np.array([ com_vx, com_vy, com_vz ])
                vel4 -= v_cen

                vz = vel4[:,2]

                vxr = xyz4 * vel4
                #v_dot_r = vxr[:,0] + vxr[:,1] + vxr[:,2]
                v_dot_r = vxr[:,0] + vxr[:,1]
                #vr = v_dot_r / ( r4 + 1.e-3 )
                vr = v_dot_r / ( r4_Cyl + 1.e-3 )

                vz_disps.append( disp(vz, m4) )
                vr_disps.append( disp(vr, m4) )

                times.append(SimTime)

            vz_disps_allruns.append(vz_disps)
            vr_disps_allruns.append(vr_disps)

    vz_disps_allruns = np.array(vz_disps_allruns)
    vr_disps_allruns = np.array(vr_disps_allruns)

    np.save('data/Fig8-VelDispZ.npy',vz_disps_allruns)
    np.save('data/Fig8-VelDispR.npy',vr_disps_allruns)

else:

    vz_disps_allruns = np.load('data/Fig8-VelDispZ.npy')
    vr_disps_allruns = np.load('data/Fig8-VelDispR.npy')

print( np.shape(vz_disps_allruns), np.shape(vr_disps_allruns) )

fig,axs = plt.subplots(1,1,figsize=(6,5.5))
plt.subplots_adjust(hspace=0.0)

print([ np.average(vz_disps_allruns[j]) for j in range(10) ], 'z velocity dispersions')
print([ np.average(vr_disps_allruns[j]) for j in range(10) ], 'r velocity dispersions')

axs.errorbar( gammas, [ np.average(vz_disps_allruns[j]) for j in range(6) ], xerr=gamma_errs, yerr=[ np.std(vz_disps_allruns[j]) for j in range(6) ], color='green', fmt='o', label='Out-of-Plane $\sigma_z$' )
axs.errorbar( gammas, [ np.average(vr_disps_allruns[j]) for j in range(6) ], xerr=gamma_errs, yerr=[ np.std(vz_disps_allruns[j]) for j in range(6) ], color='purple', fmt='o', label='In-Plane $\sigma_R$')

axs.set_xlabel('Effective Anisotropy $\gamma$', fontsize=14.5)
axs.set_ylabel('Stellar Velocity Dispersion [km/s]', fontsize=14)

xs_oh_me_oh_my = 1.0*np.ones(100)
ys_oh_me_oh_my = np.linspace(0.0,30.0,100)
axs.plot( xs_oh_me_oh_my, ys_oh_me_oh_my, color='black', linestyle='dashed', alpha=0.39 )

axs.text( 1.06, 27.3, 'Isotropic Inj.', rotation=90, va='top', color='black', alpha=0.42, fontsize=13 )
axs.text( 1.08, 28.2, r'Out-of-Disk Favored $\rightarrow$', ha='left', va='bottom', color='#553399', alpha=0.69, fontsize=13.5 )
axs.text( 0.94, 28.2, r'$\leftarrow$ In-Disk Favored', ha='right', va='bottom', color='#993355', alpha=0.69, fontsize=13.5 )

axs.tick_params(labelsize=13.5,which='both',direction='in')

axs.legend(loc='lower right',fontsize=13.5)

axs.set_xlim((0.05,20))
axs.set_ylim((0.0,30.0))

axs.xaxis.set_ticks_position('both')
axs.yaxis.set_ticks_position('both')

axs.set_xscale('log')

fig.savefig('OfficialPaperPlots/F8-VelDispByRun.png')

plt.show()
plt.close()


