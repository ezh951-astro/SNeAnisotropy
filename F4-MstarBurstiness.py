import numpy as np
import h5py
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
#from matplotlib.colors import BoundaryNorm
#from matplotlib.ticker import MaxNLocator
import zHead as hd
#from astropy.cosmology import FlatLambdaCDM

arg_options = hd.handle_arguments()
ActiveRunIDs = arg_options[0]
print(ActiveRunIDs)
date = arg_options[1]

ActiveRunIDs = [0,2,5]

###

tauavg_list = [ 5.0, 10.0, 25.0, 50.0 ]
#gamma_list = [1.0, 3.0, 10.0, 30.0, 0.3, 0.1]

gamma_dat = np.loadtxt('data/gamma_dat')

gamma_list = gamma_dat[:6,1]
gamma_errs_list = gamma_dat[:6,2]

###

#fig0,axs0 = plt.subplots(1,1,figsize=(7,5.5))

sfr_evo_array_all_runs = []
burstiness_all_runs = np.zeros((6,4))
Mstar_all_runs = np.zeros(6)
for j,(RunName, SnapDir, f) in enumerate(zip(hd.runlabels,hd.snap_dirs,hd.finalsnap_list)):

    if j in [0,1,2,3,4,5]:

        sfr_file = np.loadtxt(SnapDir[:-9]+'sfr.txt')
        time = 1000. * np.array(sfr_file[:,0][::10])
        sfr = np.array(sfr_file[:,2][::10])

        m2tot = 1.e10 * f[u'Header'].attrs[u'MassTable'][2] * len(f[u'PartType2'][u'ParticleIDs'])
        m3tot = 1.e10 * f[u'Header'].attrs[u'MassTable'][3] * len(f[u'PartType3'][u'ParticleIDs'])
        m4tot = 1.e10 * np.sum(np.array(f[u'PartType4'][u'GFM_InitialMass']))

        print( f[u'Header'].attrs[u'MassTable'] )
        print( 'm2, m3, m4', m2tot, m3tot, m4tot )
        Mstar_all_runs[j] = m4tot
        Nt = len(time)

        time_cut_array = []
        sfr_cut_array = []
        sfr_smooth_cut_array = []

        for ktau,tauavg in enumerate(tauavg_list):
            kkt = [ (time<time[i])*(time>time[i]-tauavg) for i in range(Nt) ]
            sfr_smooth = np.array([ np.average( sfr[kkt[i]] ) for i in range(Nt) ])

            kkt500 = time > 500.
            time_cut = time[kkt500]
            sfr_cut = sfr[kkt500]

            print( len(kkt500), len(sfr_smooth) )

            sfr_smooth_cut = sfr_smooth[kkt500]
            burstiness = np.std(sfr_smooth_cut) / np.average(sfr_cut)

            time_cut_array.append( time_cut )
            sfr_cut_array.append( sfr_cut )
            sfr_smooth_cut_array.append( sfr_smooth_cut )
            burstiness_all_runs[j,ktau] = burstiness

        sfr_evo_array = np.array([ time_cut_array, sfr_cut_array, sfr_smooth_cut_array ])
        sfr_evo_array_all_runs.append( sfr_evo_array )

#    if j in ActiveRunIDs:

#        axs0.plot( sfr_evo_array_all_runs[j][0,1], sfr_evo_array_all_runs[j][2,1], label=RunName )

#        axs0.set_xlabel('Time [Myr]',fontsize=14.5)
#        axs0.set_ylabel(r'SFR [${\rm M}_{\odot}$/yr]',fontsize=14.5)
#        axs0.set_xlim((500,2000))
#        axs0.set_ylim((0,3.e-1))

#        axs0.tick_params(labelsize=13)

#        axs0.legend(loc='upper right', fontsize=14)

#fig0.show()
#plt.savefig('plots/sfr_all.png')

###

fig,axs = plt.subplots(2,1,figsize=(6,6.5))
plt.subplots_adjust(hspace=0.0)

axs[0].errorbar( gamma_list, np.log10(Mstar_all_runs), xerr = gamma_errs_list, fmt='o', color='red' )

print(np.log10(Mstar_all_runs))

axs[1].errorbar( gamma_list, burstiness_all_runs[:,0], xerr = gamma_errs_list, fmt='o', color='#33bb33', alpha=0.5, label=r'$\tau_{\rm avg}$ = 5 Myr')
axs[1].errorbar( gamma_list, burstiness_all_runs[:,1], xerr = gamma_errs_list, fmt='o', color='#338833', alpha=0.5, label=r'$\tau_{\rm avg}$ = 10 Myr' )
axs[1].errorbar( gamma_list, burstiness_all_runs[:,2], xerr = gamma_errs_list, fmt='o', color='#335533', alpha=0.5, label=r'$\tau_{\rm avg}$ = 25 Myr' )

y0s = np.linspace(7.25,7.85,100)
y1s = np.linspace(0.135,1.225,100)
xs_oh_me_oh_my = 1.0 * np.ones(100)

axs[0].plot( xs_oh_me_oh_my, y0s, linestyle='dashed', alpha=0.42, color='black' )
axs[1].plot( xs_oh_me_oh_my, y1s, linestyle='dashed', alpha=0.42, color='black' )

axs[0].text( 1.05, 7.28, 'Isotropic Inj.', rotation=90, va='bottom', color='black', alpha=0.42, fontsize=12.5 )
axs[1].text( 1.07, 0.135, r'Out-of-Disk Favored $\rightarrow$', ha='left', va='bottom', color='#553399', alpha=0.69, fontsize=13 )
axs[1].text( 0.94, 0.135, r'$\leftarrow$ In-Disk Favored', ha='right', va='bottom', color='#993355', alpha=0.69, fontsize=13 )

axs[0].set_xscale('log')
axs[1].set_xscale('log')

axs[1].set_xlabel('Effective Anisotropy $\gamma$', fontsize=13.5)
axs[0].set_ylabel(r'${\rm log}_{10}$ $M_*$ Formed [${\rm M}_{\odot}$]', fontsize=13.5)
axs[1].set_ylabel(r'SFR Burstiness B', fontsize=14)

axs[0].set_xticks([])
axs[0].tick_params( labelsize=13, which='both', direction='in' )

axs[1].tick_params( labelsize=13, which='both', direction='in' )

axs[1].legend(loc='upper right')

axs[0].set_xlim((0.05,20))
axs[1].set_xlim((0.05,20))

axs[0].set_ylim((7.25,7.85))
axs[1].set_ylim((0.135,1.225))

axs[0].xaxis.set_ticks_position('both')
axs[1].xaxis.set_ticks_position('both')
axs[0].yaxis.set_ticks_position('both')
axs[1].yaxis.set_ticks_position('both')

fig.savefig('OfficialPaperPlots/F4-BurstinessMstar.png')

plt.show()
#fig.show()
