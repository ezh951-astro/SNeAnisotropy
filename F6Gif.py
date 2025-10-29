import numpy as np
import h5py

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

import TorreyLabTools.visualization.contour_makepic as makepic
import TorreyLabTools.util.calc_hsml as calc_hsml

import imageio.v2 as imageio

import zHead as hd

arg_options = hd.handle_arguments()
ActiveRunIDs = [5,0,2]
print(ActiveRunIDs)
date = arg_options[1]

###

MODE = hd.DEFAULT_MODE

RMAX = 5. # kpc
NPIX = 501
edge_aspect_ratio = 0.4

GAS_LO = -2.0
GAS_HI = 3.0

STAR_LO = -3.0
STAR_HI = 2.0

FSNAP = 130
LSNAP = 370

###

cmapGas, normGas = hd.create_cmap('magma',lo=GAS_LO,hi=GAS_HI,Nbins=50)
cmapStar, normStar = hd.create_cmap('viridis',lo=STAR_LO,hi=STAR_HI,Nbins=50)

###

Nruns = len(hd.runlabels)
ImgList = []
for snap_i in range(FSNAP,LSNAP):

      try:

            fig0,axs0 = plt.subplots(2,3,figsize=(12,5.5),gridspec_kw={ 'width_ratios':[1,1,1], 'height_ratios':[1.0,edge_aspect_ratio] })
            fig0.subplots_adjust(hspace=0.01,wspace=0.01,left=0.06,right=0.86,bottom=0.05,top=0.88)

            #fig4,axs4 = plt.subplots(2,3,figsize=(12,5.5),gridspec_kw={ 'width_ratios':[1,1,1], 'height_ratios':[1.0,edge_aspect_ratio] })
            #fig4.subplots_adjust(hspace=0.01,wspace=0.01,left=0.06,right=0.86,bottom=0.05,top=0.88)

            for k,runid in enumerate(ActiveRunIDs):

                  print(snap_i,runid)

                  RunName = hd.runlabels[runid]
                  SnapDir = hd.snap_dirs[runid]

                  print(SnapDir)

                  f = h5py.File(SnapDir + hd.get_istr3(snap_i) + '.hdf5')

                  xyz0 = np.array(f[u'PartType0'][u'Coordinates'])
                  m0 = np.array(f[u'PartType0'][u'Masses'])

                  xyz4 = np.array(f[u'PartType4'][u'Coordinates'])
                  m4 = np.array(f[u'PartType4'][u'Masses'])

                  com = hd.calc_com( xyz4, m4 )
                  r0 = hd.calc_distance( xyz0, com )
                  #r4 = hd.calc_distance( xyz4, com )
                  #rHalf = hd.half_mass_radius( r4, m4 )

                  xyz0 -= com
                  #xyz4 -= com

                  ###

                  gas_hsml = calc_hsml.get_particle_hsml( xyz0[:,0], xyz0[:,1], xyz0[:,2], DesNgb=16 )

                  massmap_0f,face_image_gas = makepic.contour_makepic( xyz0[:,0], xyz0[:,1], xyz0[:,2], gas_hsml, m0,
                        xlen = RMAX,
                        pixels = NPIX, set_aspect_ratio = 1.0,
                        set_maxden = 10.0**GAS_HI/1.e4, ## (gadget units, 10^10 msun/kpc^2 = 10^4 msun/pc^2)
                        set_dynrng = 10.0**(GAS_HI-GAS_LO)  )

                  massmap_0e,edge_image_gas = makepic.contour_makepic( xyz0[:,0], xyz0[:,2], xyz0[:,1], gas_hsml, m0,
                        xlen = RMAX,
                        pixels = NPIX, set_aspect_ratio = edge_aspect_ratio,
                        set_maxden = 10.0**GAS_HI/1.e4, ## (gadget units, 10^10 msun/kpc^2 = 10^4 msun/pc^2)
                        set_dynrng = 10.0**(GAS_HI-GAS_LO)  )

                  edge_image_gas = np.transpose(edge_image_gas) 
                  face_image_gas = np.array(face_image_gas)

                  ###

                  #stars_hsml = calc_hsml.get_particle_hsml( xyz4[:,0], xyz4[:,1], xyz4[:,2], DesNgb=32 )

                  #massmap_4f,face_image_stars = makepic.contour_makepic( xyz4[:,0], xyz4[:,1], xyz4[:,2], stars_hsml, m4,
                        #xlen = RMAX,
                        #pixels = NPIX, set_aspect_ratio = 1.0,
                        #set_maxden = 10.0**STAR_HI/1.e4 , ## (gadget units, 10^10 msun/kpc^2 = 10^4 msun/pc^2)
                        #set_dynrng = 10.0**(STAR_HI-STAR_LO)  )

                  #massmap_4e,edge_image_stars = makepic.contour_makepic( xyz4[:,0], xyz4[:,2], xyz4[:,1], stars_hsml, m4,
                        #xlen = RMAX,
                        #pixels = NPIX, set_aspect_ratio = edge_aspect_ratio,
                        #set_maxden = 10.0**STAR_HI/1.e4, ## (gadget units, 10^10 msun/kpc^2 = 10^4 msun/pc^2)
                        #set_dynrng = 10.0**(STAR_HI-STAR_LO)  )

                  #edge_image_stars = np.transpose(edge_image_stars) 
                  #face_image_stars = np.array(face_image_stars)

                  ###

                  csfac_gas = (GAS_HI-GAS_LO) / 256.
                  csfac_star = (STAR_HI-STAR_LO) / 256.

                  imFace0 = axs0[0,k].imshow( face_image_gas * csfac_gas + GAS_LO, extent=(-RMAX,RMAX, -RMAX,RMAX), origin='lower', cmap=cmapGas, norm=normGas)
                  imEdge0 = axs0[1,k].imshow( edge_image_gas * csfac_gas + GAS_LO, extent=(-RMAX,RMAX, -RMAX*edge_aspect_ratio, RMAX*edge_aspect_ratio), origin='lower', cmap=cmapGas, norm=normGas)

                  #imFace4 = axs4[0,k].imshow( face_image_stars * csfac_star + STAR_LO, extent=(-RMAX,RMAX, -RMAX,RMAX), origin='lower', cmap=cmapStar, norm=normStar)
                  #imEdge4 = axs4[1,k].imshow( edge_image_stars * csfac_star + STAR_LO, extent=(-RMAX,RMAX, -RMAX*edge_aspect_ratio, RMAX*edge_aspect_ratio), origin='lower', cmap=cmapStar, norm=normStar)

                  # Title and Ticks

                  axs0[0,0].set_title('In-Disk Favored', fontsize=18)
                  axs0[0,1].set_title('Isotropic', fontsize=18)
                  axs0[0,2].set_title('Out-of-Disk Favored', fontsize=18)

                  #axs4[0,0].set_title('In-Disk Favored', fontsize=18)
                  #axs4[0,1].set_title('Isotropic', fontsize=18)
                  #axs4[0,2].set_title('Out-of-Disk Favored', fontsize=18)

                  for l in range(2):
                        axs0[l,k].set_xticks([])
                        axs0[l,k].set_yticks([])
                        #axs4[l,k].set_xticks([])
                        #axs4[l,k].set_yticks([])

                  ### Half Mass Radius

                  #theta = np.linspace(0,2*np.pi,1000)
                  #xcirc = rHalf * np.cos(theta)
                  #ycirc = rHalf * np.sin(theta)

                  #print( 'rHalf', rHalf )

                  #for l in range(2):
                        #axs4[l,k].plot(xcirc,ycirc,linestyle='dashed',color='white')
                        #axs4[l,k].plot([0],[0], marker='+',color='white')

            # Color Bar

            cbar0_ax = fig0.add_axes([0.87,0.06,0.03,0.80])
            cbar0 = fig0.colorbar( imEdge0, cax=cbar0_ax, shrink=0.9, ticks=[-2,-1,0,1,2,3])
            cbar0.set_label(r'Gas Column Density [${\rm M}_{\odot}$ ${\rm pc}^{-2}$]', size=16)
            cbar0.ax.set_yticklabels([ '$10^{-2}$', '$10^{-1}$', '$10^0$', '$10^1$', '$10^2$', '$10^3$' ])

            #cbar4_ax = fig4.add_axes([0.87,0.06,0.03,0.80])
            #cbar4 = fig4.colorbar( imEdge4, cax=cbar4_ax, shrink=0.9, ticks=[-3,-2,-1,0,1,2])
            #cbar4.set_label(r'Stellar Column Density [${\rm M}_{\odot}$ ${\rm pc}^{-2}$]', size=16)
            #cbar4.ax.set_yticklabels([ '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^0$', '$10^1$', '$10^2$' ])

            #for cbar in [cbar0,cbar4]:
            cbar0.ax.tick_params(labelsize=15)

            # Scale Bar

            feet = fm.FontProperties(size=18)

            scalebar0 = AnchoredSizeBar( axs0[1,0].transData, 1, '1 kpc', 'lower left', pad=0.15, color='white', frameon=False, size_vertical=0.05, fontproperties=feet )
            axs0[0,0].add_artist(scalebar0)

            #scalebar4 = AnchoredSizeBar( axs4[1,0].transData, 1, '1 kpc', 'lower left', pad=0.15, color='white', frameon=False, size_vertical=0.05, fontproperties=feet )
            #axs4[0,0].add_artist(scalebar4)

            fig0.savefig('zImgTmp/imageMorph.png')
            ImgList.append( imageio.imread('zImgTmp/imageMorph.png') )

            plt.close()

      except FileNotFoundError:

            print('File Not Found')
            pass

imageio.mimsave( 'Gas.gif', ImgList )

