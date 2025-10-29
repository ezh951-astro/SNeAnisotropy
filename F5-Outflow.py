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

FSNAP = 0
LSNAP = 400

CALCULATE = False
#SFGAS_INRAD = False

YMAX_PLOT = 8.5
YMIN_PLOT = -2

###

ActiveRunIDs = [0,2,5]

print(hd.runlabels)

for j,(RunName, SnapDir) in enumerate(zip(hd.runlabels,hd.snap_dirs)):

    if CALCULATE and (j in ActiveRunIDs):

        NSNAPS = LSNAP+1

        times = np.zeros(NSNAPS)

        sfr_rad = np.zeros(NSNAPS)

        outflow_R = np.zeros(NSNAPS)
        outflow_Z = np.zeros(NSNAPS)

        outflow_MLR = np.zeros(NSNAPS)
        outflow_MLZ = np.zeros(NSNAPS)

        sfgasvel_R = np.zeros(NSNAPS)
        sfgasvel_Z = np.zeros(NSNAPS)

        sfgas_totmass = np.zeros(NSNAPS)

        for i in range(FSNAP,NSNAPS):

            print('doing {} from {} to {}'.format(i,FSNAP,LSNAP))

            f = h5py.File(SnapDir + hd.get_istr3(i) + '.hdf5')

            SimTime = 1000. * f[u'Header'].attrs[u'Time']
            times[i] = SimTime

            ##

            xyz0 = np.array(f[u'PartType0'][u'Coordinates'])
            vel0 = np.array(f[u'PartType0'][u'Velocities'])
            m0 = 1.e10 * np.array(f[u'PartType0'][u'Masses'])
            CellSFR0 = np.array(f[u'PartType0'][u'StarFormationRate'])

            try:
                xyz4 = np.array(f[u'PartType4'][u'BirthPos'])
                vel4 = np.array(f[u'PartType4'][u'Velocities'])
                m4 = 1.e10 * np.array(f[u'PartType4'][u'GFM_InitialMass'])
                age4 = SimTime - 1000. * np.array(f[u'PartType4'][u'GFM_StellarFormationTime'])
            except KeyError:
                xyz4 = np.zeros((2,3))
                vel4 = np.zeros((2,3))
                m4 = np.zeros(2)
                age4 = 2013. * np.ones(2)

            #dens = np.log10( 444.444 * np.array(f[u'PartType0'][u'Density']) )

            # center of mass crap ugh so annoying

            m2Curr = 1.e10 * np.ones_like(f[u'PartType2'][u'ParticleIDs']) * f[u'Header'].attrs[u'MassTable'][2]
            m3Curr = 1.e10 * np.ones_like(f[u'PartType3'][u'ParticleIDs']) * f[u'Header'].attrs[u'MassTable'][3]
            try:
                m4Curr = 1.e10 * np.array(f[u'PartType4'][u'Masses'])
            except KeyError:
                m4Curr = np.zeros(2)
            mcm = np.concatenate((m2Curr, m3Curr, m4Curr))

            xyz2Curr = np.array(f[u'PartType2'][u'Coordinates'])
            xyz3Curr = np.array(f[u'PartType3'][u'Coordinates'])
            try:
                xyz4Curr = np.array(f[u'PartType4'][u'Coordinates'])
            except KeyError:
                xyz4Curr = np.zeros((2,3))
            xyzcm = np.concatenate((xyz2Curr, xyz3Curr, xyz4Curr))

            vel2Curr = np.array(f[u'PartType2'][u'Velocities'])
            vel3Curr = np.array(f[u'PartType3'][u'Velocities'])
            try:
                vel4Curr = np.array(f[u'PartType4'][u'Velocities'])
            except KeyError:
                vel4Curr = np.zeros((2,3))
            velcm = np.concatenate((vel2Curr, vel3Curr, vel4Curr))

            Mstar = np.sum(mcm)
            com = np.sum( np.array([ mcm * xyzcm[:,b] for b in range(3) ]) / Mstar, axis=1 )
            vcom = np.sum( np.array([ mcm * velcm[:,b] for b in range(3) ]) / Mstar, axis=1 )

            xyz0 -= com
            xyz4 -= com
            vel0 -= vcom
            vel4 -= vcom

            # actual outflow stuff

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

            kk0_SFING = (CellSFR0 > 0.)

            kk0_INRAD = (rCyl0 < REF_RAD) * (Z0 < REF_HEIGHT)
            kk4_INRAD = (rCyl4 < REF_RAD) * (Z4 < REF_HEIGHT)
            #print(sum(kk0_SFING),len(kk0_SFING))

            kk4 = ( rCyl4 < REF_RAD ) * ( Z4 < REF_HEIGHT ) * (age4 < 5.)

            outflow_R[i] = np.sum( pr0[kkCyl0] ) * ( 525600. * 60. ) / ( THICKNESS * 3.e16 )
            outflow_Z[i] = np.sum( pz0[kkZ0] ) * ( 525600. * 60. ) / ( THICKNESS * 3.e16 )
            sfr_rad[i] = np.sum( m4[kk4] ) / ( TAU_AVG_SFR * 1.e6 ) # 1 and 2 kpc just being used as a proxy for in and out

            sfgasvel_R[i] = np.sum( pr0[kk0_SFING] ) / np.sum( m0[kk0_SFING] )
            sfgasvel_Z[i] = np.sum( pz0[kk0_SFING] ) / np.sum( m0[kk0_SFING] )

            sfgas_totmass[i] = np.sum( m0[kk0_SFING] )

            # mass loading

            kk_MLR = (pr0 > 0.) * kkCyl0
            outflow_MLR[i] = np.sum( pr0[kk_MLR] ) * ( 525600. * 60. ) / ( THICKNESS * 3.e16 )

            kk_MLZ = (pz0 > 0.) * kkZ0
            outflow_MLZ[i] = np.sum( pz0[kk_MLZ] ) * ( 525600. * 60. ) / ( THICKNESS * 3.e16 )

            cylfrac_gas = np.sum(m0[kk0_INRAD]) / np.sum(m0)
            cylfrac_stars = np.sum(m4[kk4_INRAD]) / np.sum(m4)
            cylfrac_baryons = ( np.sum(m0[kk0_INRAD]) + np.sum(m4[kk4_INRAD]) ) / ( np.sum(m0) + np.sum(m4) )
            cylfrac_SFgas = np.sum(m0[kk0_INRAD*kk0_SFING]) / np.sum(m0[kk0_SFING])

            #print('Fraction of gas in cylinder: {}'.format( cylfrac_gas ))
            if cylfrac_stars < 0.4 or cylfrac_stars > 0.6:
                print('Fraction of stars in cylinder: {}'.format( cylfrac_stars ))
            if cylfrac_baryons < 0.05 or cylfrac_baryons > 0.15:
                print('Fraction of baryons in cylinder: {}'.format( cylfrac_baryons ))
            #print('Fraction of SF gas in cylinder: {}'.format( cylfrac_SFgas ))

        mass_loading_factor_R = np.sum(outflow_MLR[outflow_MLR>0.]) / np.sum(sfr_rad)
        mass_loading_factor_Z = np.sum(outflow_MLZ[outflow_MLR>0.]) / np.sum(sfr_rad)
        print('For run {} R mass loading is {}'.format(j,mass_loading_factor_R))
        print('For run {} Z mass loading is {}'.format(j,mass_loading_factor_Z))

        OutflowRingData = np.transpose(np.array([ times, outflow_R, outflow_Z, sfr_rad, sfgasvel_R, sfgasvel_Z, sfgas_totmass ]))
        np.save('data/OutflowRingData_run{}.npy'.format(j),OutflowRingData)

fig,axs = plt.subplots(3,1,figsize=(9.5,9))
plt.subplots_adjust(hspace=0.23)

dat_z1 = np.load('data/OutflowRingData_run0.npy')
dat_z10 = np.load('data/OutflowRingData_run2.npy')
dat_z01 = np.load('data/OutflowRingData_run5.npy')

axs[0].plot( dat_z01[:,0], dat_z01[:,1], color='blue', alpha=0.69, linewidth=1 )
axs[0].plot( dat_z01[:,0], dat_z01[:,2], color='green', alpha=0.69, linewidth=1 )
ax0star = axs[0].twinx()
ax0star.plot( dat_z01[:,0], dat_z01[:,3], color='red', alpha=0.69, linewidth=1 )
ax0star.set_ylim(( YMIN_PLOT/25., YMAX_PLOT/25. ))

special_blue_1 = np.where( dat_z1[:,1] < 8.1, dat_z1[:,1], 8.15*np.ones_like(dat_z1[:,1]) )
axs[1].plot( dat_z1[:,0], special_blue_1, color='blue', alpha=0.69, linewidth=1 )
axs[1].plot( dat_z1[:,0], dat_z1[:,2], color='green', alpha=0.69, linewidth=1 )
ax1star = axs[1].twinx()
ax1star.plot( dat_z1[:,0], dat_z1[:,3], color='red', alpha=0.69, linewidth=1 )
ax1star.set_ylim(( YMIN_PLOT/25., YMAX_PLOT/25. ))

axs[2].plot( dat_z10[:,0], dat_z10[:,1], color='blue', alpha=0.69, linewidth=1, label = r'Radial Flux (${\rm r}=1$ kpc)' )
axs[2].plot( dat_z10[:,0], dat_z10[:,2], color='green', alpha=0.69, linewidth=1, label = r'Vertical Flux (${\rm z}=\pm 0.25$ kpc)' )
ax2star = axs[2].twinx()
ax2star.plot( dat_z10[:,0], dat_z10[:,3], color='red', alpha=0.69, linewidth=1, label = 'SFR' )
ax2star.set_ylim(( YMIN_PLOT/25., YMAX_PLOT/25. ))

axs[2].legend( loc = 'upper left', fontsize=15.5 )
ax2star.legend( loc = 'upper right', fontsize=15.5 )

axs[0].set_xticklabels([])
axs[1].set_xticklabels([])

for k in range(3):
    axs[k].set_xlim((0.,2000.))
    axs[k].set_ylim((YMIN_PLOT,YMAX_PLOT))
    axs[k].xaxis.set_ticks_position('both')
    axs[k].yaxis.set_ticks_position('left')
    axs[k].set_yticks([0,3,6])
    axs[k].tick_params(labelsize=15.5,direction='in')

for axstar in [ax0star,ax1star,ax2star]:
    axstar.set_yticks([0,0.3])
    axstar.tick_params(labelsize=15.5,direction='in')

#ax0star.tick_params(labelsize=14,direction='in')
#ax1star.tick_params(labelsize=14,direction='in')
#ax2star.tick_params(labelsize=14,direction='in')

axs[0].set_title('In-Disk Favored',fontsize=18)
axs[1].set_title('Isotropic',fontsize=18)
axs[2].set_title('Out-of-Disk Favored',fontsize=18)

#fig.supxlabel('Time [Myr]',fontsize=18)
#fig.supylabel(r'Mass Flux [${\rm M}_{\odot}$/yr]',fontsize=18)

ax1star.set_ylabel(r'SFR within Cylinder [${\rm M}_{\odot}$ ${\rm yr}^{-1}$]',fontsize=16.5)
axs[1].set_ylabel(r'Mass Flux [${\rm M}_{\odot}$ ${\rm yr}^{-1}$]',fontsize=16.5)
axs[2].set_xlabel('Time [Myr]', fontsize=16.5)

#fig.text( 0.93, 0.51, r'SFR [${\rm M}_{\odot}/yr$]', rotation=90, fontsize=16, ha='center', va='center')
#fig.text( 0.03, 0.51, r'Mass Flux [${\rm M}_{\odot}/yr$]', rotation=90, fontsize=16, ha='center', va='center')
#fig.text( 0.48, 0.05, 'Time [Myr]', fontsize=16, ha='center')

fig.savefig('OfficialPaperPlots/F5-OutflowHistory.png')
plt.show()
