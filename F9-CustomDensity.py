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

MODE = hd.DEFAULT_MODE
ActiveRunIDs = [0,2,5]
COLORS = ['purple',1,'C1',3,4,'green']

Nrad = 500

rad_log = 10.**np.linspace(-1.5,1.5,Nrad)

FSNAP = 200
LSNAP = 400
N_INIT = 40

NSNAP = LSNAP - FSNAP + 1

DYN_DIST = 2.0 # kpc

CALCULATE = True

fig,axs = plt.subplots(1,1,figsize=(8.75,6.25))

if CALCULATE:

    for j,(RunName, SnapDir) in enumerate(zip(hd.runlabels,hd.snap_dirs)):
        if j in ActiveRunIDs:

            print(j)

            all_densities = []
            for i in range(401):

                f = h5py.File(SnapDir + hd.get_istr3(i) + '.hdf5')

                m2 = 1.e10 * f[u'Header'].attrs[u'MassTable'][2] * np.ones_like(np.array(f[u'PartType2'][u'ParticleIDs']))
                m3 = 1.e10 * f[u'Header'].attrs[u'MassTable'][3] * np.ones_like(np.array(f[u'PartType3'][u'ParticleIDs']))

                try:
                    m4 = 1.e10 * np.array(f[u'PartType4'][u'Masses'])
                except KeyError:
                    m4 = np.zeros(2)
                mAS = np.concatenate((m2,m3,m4))

                xyz2 = np.array(f[u'PartType2'][u'Coordinates'])
                xyz3 = np.array(f[u'PartType3'][u'Coordinates'])
                try:
                    xyz4 = np.array(f[u'PartType4'][u'Coordinates'])
                except KeyError:
                    xyz4 = np.zeros((2,3))
                xyzAS = np.concatenate((xyz2,xyz3,xyz4))

                Mstar = np.sum(mAS)
                com = np.sum( np.array([ mAS*xyzAS[:,b] for b in range(3) ])/Mstar, axis=1 )
                print( j,i,com )

                MDlog = hd.mass_distribution_spherical(f, rad_log, mode=MODE, origin=[ com[b] for b in range(3) ], ptypes=[0,1] )

                rads = MDlog[1,0]
                density = MDlog[1,2]

                # For tdyn:

                argdyn = np.argmax( rads > DYN_DIST )
                raddyn = rads[ argdyn ]
                densdyn = MDlog[0,2][argdyn] + MDlog[1,2][argdyn]
                tdyn = np.sqrt( 3.0 / 4.0 / np.pi / hd.G / densdyn ) * 951. # num of Myrs in a kpc/(km/s)
                print( 'Dynamical Timescale At {} kpc is {} Myr'.format(raddyn,tdyn) )

                all_densities.append(density)
                #print(density[::25])

            all_densities = np.array(all_densities)
            print(np.shape(all_densities))
            np.savetxt('data/alldensity{}.txt'.format(j),all_densities)

all_densities = np.loadtxt('data/alldensity5.txt')
# initial density
median_initial_density = []
for krad in range(Nrad):
    init_dens_ranked = np.sort(all_densities[:N_INIT+1,krad])
    median_initial_density.append( init_dens_ranked[ int( 0.5*(N_INIT+1) ) ] )

axs.plot( rad_log, np.array(median_initial_density)/1.e9, label='Initial Density', linewidth=2.0, color='#999999', alpha=0.69, linestyle='dashed' )

for j,(RunName, SnapDir) in enumerate(zip(hd.runlabels,hd.snap_dirs)):
    if j in ActiveRunIDs:

        all_densities = np.loadtxt('data/alldensity{}.txt'.format(j))

        percentiles = np.zeros((Nrad,4))
        for krad in range(Nrad):
            dens_ranked = np.sort(all_densities[FSNAP:LSNAP+1,krad])
            p25th = dens_ranked[ int( 0.25*NSNAP ) ]
            p50th = dens_ranked[ int( 0.50*NSNAP ) ]
            p75th = dens_ranked[ int( 0.75*NSNAP ) ]
            percentiles[krad,1] = p25th
            percentiles[krad,2] = p50th
            percentiles[krad,3] = p75th

        percentiles[:,0] = rad_log

        np.savetxt('data/densitypercent{}.txt'.format(j),percentiles)

        axs.fill_between( percentiles[:,0], percentiles[:,1]/1.e9, y2=percentiles[:,3]/1.e9, alpha=0.2, color=COLORS[j] )
        axs.plot( percentiles[:,0], percentiles[:,2]/1.e9, label=RunName, linewidth=2.5, color=COLORS[j], alpha=0.54 )

        # Dynamical Timescale of the Inner Kiloparsec

        #DYN_DIST = 1.0

        #dens25 = percentiles[ argdyn, 1 ]
        #dens75 = percentiles[ argdyn, 3 ]

        #tdyn25 = np.sqrt( 3.0 / 4.0 / np.pi / hd.G / dens25 ) * 951. # num of Myrs in a kpc/(km/s)
        #tdyn75 = np.sqrt( 3.0 / 4.0 / np.pi / hd.G / dens75 ) * 951. # num of Myrs in a kpc/(km/s)

        #print( '25th: tdyn={} Myr'.format(tdyn25) )
        #print( '75th: tdyn={} Myr'.format(tdyn75) )

        axs.legend(loc='upper right',fontsize=17)
        axs.set_xlabel('Radius [kpc]',fontsize=17.5)
        axs.set_ylabel('Dark Density [$\mathrm{M}_{\odot}$ $\mathrm{pc^{-3}}$]',fontsize=17.5)

        axs.set_xscale('log')
        axs.set_yscale('log')

        axs.set_xlim((0.1,3.0))
        axs.set_xticks([0.1, 0.3, 1.0, 3.0])
        axs.set_xticklabels([0.1, 0.3, 1.0, 3.0])

        axs.set_ylim((0.02,1.0))
        axs.set_yticks([0.03, 0.1, 0.3, 1.0])
        axs.set_yticklabels([0.03, 0.1, 0.3, 1.0])

        axs.xaxis.set_ticks_position( 'both' )
        axs.yaxis.set_ticks_position( 'both' )

        axs.tick_params(labelsize=16.5, which='both', direction='in')

fig.savefig('OfficialPaperPlots/F9-DMDensity.png')
plt.show()
