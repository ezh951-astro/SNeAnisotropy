import numpy as np
import h5py
import matplotlib.pyplot as plt
import zHead as hd
from astropy.cosmology import FlatLambdaCDM

arg_options = hd.handle_arguments()
ActiveRunIDs = arg_options[0]
print(ActiveRunIDs)
date = arg_options[1]
NRUNS = len(ActiveRunIDs)

###

MODE = hd.DEFAULT_MODE

REF_RAD = 1.0
REF_HEIGHT = 0.25

THICKNESS = 0.05

TAU_AVG_SFR = 5. # Myr; fixed

FSNAP = 30
LSNAP = 300

###

CALCULATE = True

def linf(x,m,b):
    return m*x+b

###

print(hd.runlabels)

if CALCULATE:
    for j,(RunName, SnapDir) in enumerate(zip(hd.runlabels,hd.snap_dirs)):
        if j in ActiveRunIDs:

            NSNAPS = LSNAP+1
            times = np.zeros(NSNAPS)
            sfr_rad = np.zeros(NSNAPS)
            outflow_R = np.zeros(NSNAPS)
            outflow_Z = np.zeros(NSNAPS)

            for i in range(FSNAP,NSNAPS):

                print('doing {} from {} to {}'.format(i,FSNAP,LSNAP))

                f = h5py.File(SnapDir + hd.get_istr3(i) + '.hdf5')

                SimTime = 1000. * f[u'Header'].attrs[u'Time']
                times[i] = SimTime

                ##

                xyz0 = np.array(f[u'PartType0'][u'Coordinates'])
                xyz4 = np.array(f[u'PartType4'][u'BirthPos'])

                vel0 = np.array(f[u'PartType0'][u'Velocities'])
                vel4 = np.array(f[u'PartType4'][u'Velocities'])

                m0 = 1.e10 * np.array(f[u'PartType0'][u'Masses'])
                m4 = 1.e10 * np.array(f[u'PartType4'][u'Masses'])

                dens = np.log10( 444.444 * np.array(f[u'PartType0'][u'Density']) )

                age4 = SimTime - 1000. * np.array(f[u'PartType4'][u'GFM_StellarFormationTime'])

                Mstar = np.sum(m4)
                com = np.sum( np.array([ m4 * xyz4[:,b] for b in range(3) ]) / Mstar, axis=1 )
                vcom = np.sum( np.array([ m4 * vel4[:,b] for b in range(3) ]) / Mstar, axis=1 )

                xyz0 -= com
                xyz4 -= com
                vel0 -= vcom
                vel4 -= vcom

                rCyl0 = np.sqrt( xyz0[:,0]*xyz0[:,0] + xyz0[:,1]*xyz0[:,1] )
                rCyl4 = np.sqrt( xyz4[:,0]*xyz4[:,0] + xyz4[:,1]*xyz4[:,1] )

                Z0 = np.abs( xyz0[:,2] )
                Z4 = np.abs( xyz4[:,2] )

                vr0 = ( xyz0[:,0]*vel0[:,0] + xyz0[:,1]*vel0[:,1] ) / (rCyl0 + 1.e-5)

                vz0plus = xyz0[:,2]*vel0[:,2] / (xyz0[:,2] + 1.e-5) # note; this is the component of the velocity in the POSITIVE Z direction.
                vz0minus = -xyz0[:,2]*vel0[:,2] / (xyz0[:,2] + 1.e-5) # note; this is the component of the velocity in the NEGATIVE Z direction.
                vz0 = np.where( xyz0[:,2] > 0., vz0plus, vz0minus )

                pr0 = m0 * vr0
                pz0 = m0 * vz0

                kkCyl0_in = rCyl0 > REF_RAD - (THICKNESS/2.)
                kkCyl0_out = rCyl0 < REF_RAD + (THICKNESS/2.)
                kkCyl0 = kkCyl0_in * kkCyl0_out * ( Z0 < REF_HEIGHT )

                kkZ0_in = Z0 > REF_HEIGHT - (THICKNESS/2.)
                kkZ0_out = Z0 < REF_HEIGHT + (THICKNESS/2.)
                kkZ0 = kkZ0_in * kkZ0_out * ( rCyl0 < REF_RAD )

                kk4 = ( rCyl4 < REF_RAD ) * ( Z4 < REF_HEIGHT )

                outflow_R[i] = np.sum( pr0[kkCyl0] ) * ( 525600. * 60. ) / ( THICKNESS * 3.e16 )
                outflow_Z[i] = np.sum( pz0[kkZ0] ) * ( 525600. * 60. ) / ( THICKNESS * 3.e16 )
                sfr_rad[i] = np.sum( m4[kk4] ) / ( TAU_AVG_SFR * 1.e6 ) # 1 and 2 kpc just being used as a proxy for in and out

            sfr_rad = sfr_rad - np.roll( sfr_rad, 1 )

            OutflowRingData = np.transpose(np.array([ times, outflow_R, outflow_Z, sfr_rad ]))
            np.save('data/OutflowRingData_run{}.npy'.format(j),OutflowRingData)

###

Params = ['HiEff','LoThr','HiEff-LoThr','HighRes']
Zetas = ['zeta1 (Isotropic)','zeta10 (Out-of-Disk Favored)']

fig,axs = plt.subplots(4,2,figsize=(9.5,8))
plt.subplots_adjust(hspace=0.0,wspace=0.05)

for j in range(8):

    a = np.load('data/OutflowRingData_run{}.npy'.format(j))

    axs[ j//2, j%2 ].plot( a[:,0], 25.*a[:,3], label='SFR x 25', color='red', alpha=0.5 )
    axs[ j//2, j%2 ].plot( a[:,0], a[:,1], label='In-Disk Outflow Rate', color='blue', alpha=0.5 )
    axs[ j//2, j%2 ].plot( a[:,0], a[:,2], label='Out-of-Disk Outflow Rate', color='green', alpha=0.52 )

    axs[ j//2, j%2 ].set_xlim((500,1500))
    axs[ j//2, j%2 ].set_ylim((-1.75,4.25))
    axs[ j//2, j%2 ].set_xticks([600,800,1000,1200,1400])
    axs[ j//2, j%2 ].set_yticks([0,2,4])
    axs[ j//2, j%2 ].set_xticklabels([])
    axs[ j//2, j%2 ].set_yticklabels([])

    axs[ j//2, j%2 ].tick_params(which='both',direction='in')
    axs[ j//2, j%2 ].xaxis.set_ticks_position('both')
    axs[ j//2, j%2 ].yaxis.set_ticks_position('both')

for k in range(2):
    axs[3,k].set_xticklabels([600,800,1000,1200,1400],fontsize=12)
    axs[3,k].set_xlabel(Zetas[k],fontsize=13)

for k in range(4):
    axs[k,0].set_yticklabels([0,2,4],fontsize=12)
    axs[k,0].set_ylabel(Params[k],fontsize=13)

axs[0,1].legend(loc='upper right')
fig.supxlabel('Time [Myr]',fontsize=15)
fig.supylabel(r'Mass Rate [${\rm M}_{\odot}$/yr]',fontsize=15)

fig.savefig('../OfficialPaperPlots/F12-Robustness.png')
plt.show()
