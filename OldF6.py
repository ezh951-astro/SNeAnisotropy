import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import convolve2d
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import zHead as hd
from astropy.cosmology import FlatLambdaCDM
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

arg_options = hd.handle_arguments()
ActiveRunIDs = [5,0,2]
print(ActiveRunIDs)
date = arg_options[1]

###

MODE = hd.DEFAULT_MODE

CENTER = 100.
RMAX = 4.

NPIX = 200
PIX_SMEAR = 1

i = 240

###

cmapGas, normGas = hd.create_cmap('magma',lo=3,hi=9,Nbins=24)

def gauss2d(x,y,r,xcen,ycen):
    return np.exp( -0.5 * ( (x-xcen)**2. + (y-ycen)**2. ) / r**2. ) / r / r / (2.*np.pi)

def gaussian_smear( arr2d , rad = 2. ):
    gaub_arr = np.array([ [gauss2d(i,j,rad,(NPIX-1.)/2.,(NPIX-1.)/2.) for j in range(NPIX)] for i in range(NPIX) ])
    renormalization_factor = np.sum(gaub_arr)
    print(renormalization_factor)
    gaub_arr = gaub_arr / renormalization_factor
    ans = convolve2d( arr2d, gaub_arr )
    return ans

Nruns = len(hd.runlabels)

fig0,axs0 = plt.subplots(2,3,figsize=(11,8))
fig0.subplots_adjust(hspace=0.01,wspace=0.01)

for k,runid in enumerate(ActiveRunIDs):

    RunName = hd.runlabels[runid]
    SnapDir = hd.snap_dirs[runid]
    
    f = h5py.File(SnapDir + hd.get_istr3(i) + '.hdf5')

    xyz4 = np.array(f[u'PartType4'][u'Coordinates'])
    mList = np.array(f[u'PartType4'][u'Masses'])

    totmass = np.sum(mList)
    com_x = np.sum(mList * xyz4[:,0]) / totmass
    com_y = np.sum(mList * xyz4[:,1]) / totmass
    com_z = np.sum(mList * xyz4[:,2]) / totmass
    com = np.array([ com_x, com_y, com_z ])

    #xyz4_com = xyz4 - com
    #rList = np.sqrt(sum([ xyz4_com[:,b]*xyz4_com[:,b] for b in range(3) ]))
    #halfMassRad4 = hd.half_mass_radius( rList, mList )

    SimTime = f[u'Header'].attrs[u'Time']
    morphology_i = hd.morphology(f, mode=MODE, origin=com, rad=RMAX, nbins=NPIX, ptypes=[0,4])

    if PIX_SMEAR < 0.01:
        Face0 = np.log10(morphology_i[0][0]+1.)
        Edge0 = np.log10(morphology_i[0][1]+1.)

    else:
        Face0 = gaussian_smear( np.log10(morphology_i[0][0]+1.), rad=PIX_SMEAR )[ int(NPIX/2.):int(3*NPIX/2.), int(NPIX/2.):int(3*NPIX/2.) ]
        Edge0 = gaussian_smear( np.log10(morphology_i[0][1]+1.), rad=PIX_SMEAR )[ int(NPIX/2.):int(3*NPIX/2.), int(NPIX/2.):int(3*NPIX/2.) ]

    imFace0 = axs0[0,k].imshow( Face0, extent=(-RMAX,RMAX,-RMAX,RMAX), origin='lower', cmap=cmapGas, norm=normGas )
    imEdge0 = axs0[1,k].imshow( Edge0, extent=(-RMAX,RMAX,-RMAX,RMAX), origin='lower', cmap=cmapGas, norm=normGas )

    axs0[0,0].set_ylabel('Gas: Face-On', fontsize=18)
    axs0[1,0].set_ylabel('Gas: Edge-On', fontsize=18)

    axs0[0,0].set_title('In-Disk\nFavored', fontsize=18)
    axs0[0,1].set_title('Isotropic', fontsize=18)
    axs0[0,2].set_title('Out-of-Disk\nFavored', fontsize=18)

    for l in range(2):
        axs0[l,k].set_xticks([])
        axs0[l,k].set_yticks([])

    #theta = np.linspace(0,2*np.pi,1000)
    #xcirc = halfMassRad4 * np.cos(theta)
    #ycirc = halfMassRad4 * np.sin(theta)

    #print( halfMassRad4 )

    #for l in range(2):
        #axs4[l,k].plot(xcirc,ycirc,linestyle='dashed',color='white')
        #axs4[l,k].plot([0],[0], marker='+',color='white')

    ###

feet = fm.FontProperties(size=18)

scalebar0 = AnchoredSizeBar( axs0[1,0].transData, 1, '1 kpc', 'lower left', pad=0.15, color='white', frameon=False, size_vertical=0.05, fontproperties=feet )
axs0[1,0].add_artist(scalebar0)

fig0.savefig('PaperPlots/ProjectionCustom_snap{}_gas.png'.format(i))

plt.show()
plt.close()

