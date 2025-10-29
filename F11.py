import numpy as np
import h5py
import matplotlib.pyplot as plt
import zHead as hd
import pickle

arg_options = hd.handle_arguments()
ActiveRunIDs = arg_options[0]
print(ActiveRunIDs)
date = arg_options[1]

###

colors = ['#AD1519','#004494','#006633']
CustomRunLabels = ['0','1','Solid-Angle Weight (Default)','3','4','5','6','Volume Weight','8','Mass Weight']

ColumnSet = [ [0,6,8], [2,7,9] ]

tauavg = 10.
FIRST_SNAPSHOT = 0

CALCULATE = False

###

if CALCULATE:

    all_10runs_star_data = [] # use all_10runs_star_data[ j_of_run ][ k_of_data ]
    for j,(RunName, SnapDir, f) in enumerate(zip(hd.runlabels,hd.snap_dirs,hd.finalsnap_list)):

        # Calculate total stellar mass formed

        times_Mstar = []
        totmasses_Mstar = []
        for i in range(FIRST_SNAPSHOT,hd.final_snap_nums[j]):
            print(i)

            f = h5py.File(SnapDir + hd.get_istr3(i) + '.hdf5')

            try:
                masses4 = 1.e10 * np.array(f[u'PartType4'][u'GFM_InitialMass'])
                #masses3 = f[u'Header'].attrs[u'MassTable'][3] * np.ones_like( np.array(f[u'PartType3'][u'ParticleIDs']) )
                #masses2 = f[u'Header'].attrs[u'MassTable'][2] * np.ones_like( np.array(f[u'PartType2'][u'ParticleIDs']) )
            except KeyError:
                masses4 = np.array([0,0])

            #print( masses3, masses2 )
            #TotMass = np.sum(masses4) + np.sum(masses3) + np.sum(masses2)
            TotMass = np.sum(masses4)

            times_Mstar.append( f[u'Header'].attrs[u'Time'] * 1000. )
            totmasses_Mstar.append(TotMass)

        # Calculate SFR, using a 10 Myr rolling average

        sfr_file = np.loadtxt(SnapDir[:-9]+'sfr.txt')
        time = 1000. * np.array(sfr_file[:,0][::10])
        sfr = np.array(sfr_file[:,2][::10])

        Nt = len(time)

        kkt = [ (time<time[i])*(time>time[i]-tauavg) for i in range(Nt) ]
        sfr_smooth = np.array([ np.average( sfr[kkt[i]] ) for i in range(Nt) ])

        #kkt500 = time > 0.

        #time_cut_SFR = time[kkt500]
        #sfr_smooth_cut_SFR = sfr_smooth[kkt500]

        # store all of the data

        run_star_data = [ time, sfr_smooth, times_Mstar, totmasses_Mstar ]
        all_10runs_star_data.append( run_star_data )

    with open("data/StarSizeData", "wb") as fp:
        pickle.dump(all_10runs_star_data, fp)

with open("data/StarSizeData", "rb") as fp:
    all_10runs_star_data = pickle.load(fp)

fig0,axs0 = plt.subplots(2,2,figsize=(14.5,7))
plt.subplots_adjust(hspace=0.0,wspace=0.15)

for krow in range(2):
    for kcol in range(2):

        ActiveRunIDs = ColumnSet[kcol]
        print(ActiveRunIDs)

        for kline in range(3):

            j_run = ActiveRunIDs[kline]
            time_plotline = all_10runs_star_data[ j_run ][ 2*krow ]
            value_plotline = all_10runs_star_data[ j_run ][ 2*krow+1 ]

            axs0[krow,kcol].plot( time_plotline, value_plotline, label=CustomRunLabels[j_run], color=colors[kline], alpha=0.69 )
            print(j_run)
            #print(time_plotline)
            print(np.log10(value_plotline[-1]))

        axs0[0,kcol].set_xticklabels([])
        axs0[1,1].legend(loc='lower right', fontsize=13)

        axs0[0,kcol].set_ylim((0,0.2))
        axs0[0,kcol].set_yticks([ 0, 0.05, 0.10, 0.15, 0.20 ])

        axs0[1,kcol].set_ylim((2.e5,2.e8))
        axs0[1,kcol].set_yticks([ 1.e6, 1.e7, 1.e8 ])
        axs0[1,kcol].set_yscale('log')

        axs0[krow,kcol].set_xlim((0,2000))

        axs0[krow,kcol].tick_params(labelsize=13.5)

        axs0[krow,kcol].xaxis.set_ticks_position('both')
        axs0[krow,kcol].yaxis.set_ticks_position('both')
        axs0[krow,kcol].tick_params(which='both',direction='in')

axs0[0,0].set_title(r'Isotropic', fontsize=16)
axs0[0,1].set_title(r'Out-of-Disk Favored', fontsize=16)

axs0[0,0].set_ylabel(r'SFR [${\rm M}_{\odot}$/yr]', fontsize=14)
axs0[1,0].set_ylabel(r'$M_{*}$ Formed [${\rm M}_{\odot}$]', fontsize=14)

fig0.supxlabel('Time [Myr]', fontsize=15)

fig0.savefig('OfficialPaperPlots/F11.png')
plt.show()
