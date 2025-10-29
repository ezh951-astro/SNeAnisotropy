import numpy as np
import matplotlib.pyplot as plt
import zHead as hd
import h5py

arg_options = hd.handle_arguments()
#ActiveRunIDs = arg_options[0]
ActiveRunIDs = [5,0,2]
print(ActiveRunIDs)
date = arg_options[1]

###

CustomRunLabels = ['Isotropic','1','Out-of-Disk Favored','3','4','In-Disk Favored']
CustomColors = ['purple','black','C1','black','black','green']

avj = 132
CALCULATE = True

fig,axs = plt.subplots(4,3,figsize=(13,6),gridspec_kw={ 'width_ratios': [1,1,1], 'height_ratios': [1.3,0.6,0.2,1], 'wspace':0.03, 'hspace':0 })

plt.subplots_adjust(hspace=0.0)

for kcol,j in enumerate(ActiveRunIDs):

    RunName = CustomRunLabels[j]
    SnapDir = hd.snap_dirs[j]
    Color = CustomColors[j]

    sfr_file = np.loadtxt(SnapDir[:-9]+'sfr.txt')

    time = 1000. * np.array(sfr_file[:,0])

    sfr_u = np.array(sfr_file[:,2])
    print(j,max(sfr_u))
    sfr_u_a = np.concatenate( (np.zeros(avj), sfr_u) )

    sfr = np.array([ np.average(sfr_u_a[ ind:ind+avj ]) for ind in range(len(time)) ])
    sfr = np.where( sfr<0.295, sfr, 0.295*np.ones_like(sfr) )

    if CALCULATE:
        m4tots = []
        time_by_snap = []
        for i in range(401):
            print(i)
            f = h5py.File( SnapDir + hd.get_istr3(i) + '.hdf5' )
            time_by_snap.append( f[u'Header'].attrs[u'Time'] * 1000 )
            try:
                m4 = np.array(f[u'PartType4'][u'Masses']) * 1.e10
                m4tot = np.sum(m4)
                print( 'average particle mass {}'.format( np.average(m4) ) )
                print( 'std particle mass {}'.format( np.std(m4) ) )
            except KeyError:
                m4tot = 1.e0
            m4tots.append(m4tot)
        np.savetxt('data/mstardat_run{}.txt'.format(j), np.transpose(np.array([time_by_snap,m4tots])) )
    else:
        mstardat = np.loadtxt('data/mstardat_run{}.txt'.format(j))
        time_by_snap = mstardat[:,0]
        m4tots = mstardat[:,1]

    axs[0,kcol].plot( time, sfr, label=RunName, alpha=0.84, color=Color, linewidth=1.0 )
    axs[1,kcol].plot( time, np.log10(sfr), label=RunName, alpha=0.84, color=Color, linewidth=1.0 )
    axs[3,kcol].plot( time_by_snap, np.log10(m4tots), label=RunName, alpha=0.84, color=Color, linewidth=1.0 )

    axs[0,kcol].set_title(RunName, fontsize=16)

fig.text( 0.0675, 0.64, r'SFR [${\rm M}_{\odot}$ ${\rm yr}^{-1}$]',fontsize=16,rotation=90,va='center',ha='center' )
#fig.text( 0.065, 0.76, r'SFR [${\rm M}_{\odot}$/yr]',fontsize=16,rotation=90,va='center' )
#fig.text( 0.065, 0.50, r'log10 SFR [${\rm M}_{\odot}$/yr]',fontsize=16,rotation=90,va='center')
fig.text( 0.0675, 0.21, 'Stellar Mass\n' + r'Formed [${\rm M}_{\odot}$]',fontsize=16,rotation=90,va='center',ha='center')

for kcol in range(3):

    axs[2,kcol].set_visible(False)

    #axs[3,kcol].set_yscale('log')
    axs[3,1].set_xlabel('Time [Myr]',fontsize=16)

    for krow in range(4):
        axs[krow,kcol].set_xlim((0,2000))
        axs[krow,kcol].set_xticks([500,1000,1500,2000],fontsize=14.5)
        axs[krow,kcol].tick_params(which='both',direction='in',labelsize=15)
        axs[krow,kcol].xaxis.set_ticks_position('both')
        axs[krow,kcol].yaxis.set_ticks_position('both')

    axs[0,kcol].set_ylim((0,0.3))
    axs[1,kcol].set_ylim((-2.9,-0.1))
    axs[3,kcol].set_ylim((5.25,7.75))

    axs[0,kcol].set_yticks( [ 0.1, 0.2, 0.3 ], fontsize=15)
    axs[1,kcol].set_yticks( [ -2.0, -1.0 ], fontsize=15)
    axs[3,kcol].set_yticks( [ 6.0, 7.0 ], fontsize=15)

    axs[0,0].set_yticklabels( [ 0.1, 0.2, 0.3 ], fontsize=15)
    axs[1,0].set_yticklabels( [ r'$10^{-2}$', r'$10^{-1}$' ], fontsize=15)
    axs[3,0].set_yticklabels( [ r'$10^6$', r'$10^7$' ], fontsize=15)

    for krow in range(2):
        axs[krow,kcol].set_xticklabels([])

for krow in range(4):
    for kcol in range(2):
        axs[krow,kcol+1].set_yticklabels([])

axs[3,0].set_xticks([0,500,1000,1500,2000])

    #axs[2,kcol].legend(loc='lower right', fontsize=16)

plt.savefig('OfficialPaperPlots/F3-SFR.png')

plt.show()
