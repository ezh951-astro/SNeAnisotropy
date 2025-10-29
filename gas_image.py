import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib import rc
#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#rc('text', usetex=True)

import numpy as np

import TorreyLabTools.visualization.contour_makepic as makepic
import TorreyLabTools.util.calc_hsml as calc_hsml
import TorreyLabTools.util.naming as naming

import TorreyLabTools.simread.readsnapHDF5 as rs

import os
import h5py

import zHead as hd

edge_aspect_ratio = 0.25
left_text_buffer = 0.0
top_text_buffer  = 0.125

image_size = 1
n_pixels = 501

IMG_BOX_WIDTH = 12.0 # kpc


files = [ hd.snap_dirs[runid]+'240.hdf5' for runid in [5,0,2] ]

#files = [
    #'/orange/paul.torrey/paul.torrey/Share/ezhang/dwarf_images/data/zeta0.1_1200Myr.hdf5', 
    #'/orange/paul.torrey/paul.torrey/Share/ezhang/dwarf_images/data/zeta1.0_1200Myr.hdf5', 
    #'/orange/paul.torrey/paul.torrey/Share/ezhang/dwarf_images/data/zeta10.0_1200Myr.hdf5' ]

top_text_array = [ r'In-Disk Favored', r'Isotropic', r'Out-of-Disk Favored']

header = rs.snapshot_header( files[0] ) 
boxsize = header.boxsize
num_files = header.filenum

all_im_fig = plt.figure(figsize=(image_size*len(files) + left_text_buffer, image_size+image_size*edge_aspect_ratio + top_text_buffer))

for file_index,this_file in enumerate(files):
    center = [boxsize/2.0,boxsize/2.0,boxsize/2.0]  
    gas_pos  = rs.read_block( this_file, 'POS ', parttype=0)
    gas_mass = rs.read_block( this_file, 'MASS', parttype=0)
    gas_vel  = rs.read_block( this_file, 'VEL ', parttype=0)

    for ijk in range(3):
      gas_pos[:,ijk] -= center[ijk]
      gas_pos[ gas_pos[:,ijk] > boxsize/2.0 , ijk ] -= boxsize
      gas_pos[ gas_pos[:,ijk] < -1.0*boxsize/2.0 , ijk ] += boxsize

    gas_hsml = calc_hsml.get_particle_hsml( gas_pos[:,0], gas_pos[:,1], gas_pos[:,2], DesNgb=16  )

    massmap,face_image = makepic.contour_makepic( gas_pos[:,0], gas_pos[:,1], gas_pos[:,2], gas_hsml, gas_mass ,
#          xlen = boxsize/40.0,
          xlen = IMG_BOX_WIDTH/2.0,
          pixels = n_pixels, set_aspect_ratio = 1.0,
          set_maxden = 1.0e-1, ## (gadget units, 10^10 msun/kpc^2 = 10^4 msun/pc^2)
          set_dynrng = 1.0e4  )


    massmap,edge_image = makepic.contour_makepic( gas_pos[:,0], gas_pos[:,2], gas_pos[:,1], gas_hsml, gas_mass ,
#          xlen = boxsize/40.0,
          xlen = IMG_BOX_WIDTH/2.0,
          pixels = n_pixels, set_aspect_ratio = edge_aspect_ratio,
          set_maxden = 1.0e-1, ## (gadget units, 10^10 msun/kpc^2 = 10^4 msun/pc^2)
          set_dynrng = 1.0e4  )

    edge_image = np.transpose( edge_image ) 

    face_image = np.array( face_image ) 

#    for colormap in ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
#                      'pink', 'afmhot', 'gist_heat', 'copper', 'twilight', 'twilight_shifted', 
#                      'gist_stern', 'gnuplot', 'gnuplot2', 'CMRmap',
#a                      'cubehelix', 'rainbow' ]:

    for colormap in ['magma']: 
      width = 1.0 /  (left_text_buffer + len(files)) 
      left  = left_text_buffer / (left_text_buffer + len(files))  + width*file_index

      bottom_height = edge_aspect_ratio / (top_text_buffer + 1.0 + edge_aspect_ratio)
      top_height    = 1.0 / (top_text_buffer + 1.0 + edge_aspect_ratio)
      
      bottom_bottom = 0.0
      top_bottom = bottom_height
      

      print(np.shape(face_image))
      print(np.shape(edge_image))

      im_ax = all_im_fig.add_axes([left, top_bottom , width, top_height ])
      im_ax.imshow( face_image, cmap='magma', extent=(-IMG_BOX_WIDTH/2.0, IMG_BOX_WIDTH/2.0, -IMG_BOX_WIDTH/2.0, IMG_BOX_WIDTH/2.0) )
      im_ax.xaxis.set_ticks_position('none') 
      im_ax.yaxis.set_ticks_position('none') 
      im_ax.tick_params(  axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        left=False,
        right=False,
        labelbottom=False,
        labeltop=False,
        labelleft=False,
        labelright=False )

      #im_ax.text( 0.0, boxsize/40.0, r'$\mathrm{Test}$' , horizontalalignment='center')
      im_ax.text( 0.0, boxsize/40.0*1.03, top_text_array[file_index] , horizontalalignment='center', fontsize=7.0)

      im_ax.spines['bottom'].set_color('white')
      im_ax.spines['top'].set_color('white')
      im_ax.spines['left'].set_color('white')
      im_ax.spines['right'].set_color('white')

      im_ax = all_im_fig.add_axes([left, bottom_bottom , width, bottom_height])
      im_ax.imshow( edge_image, cmap='magma', extent=(-IMG_BOX_WIDTH/2.0, IMG_BOX_WIDTH/2.0, -IMG_BOX_WIDTH*edge_aspect_ratio/2.0, IMG_BOX_WIDTH*edge_aspect_ratio/2.0) )

      im_ax.xaxis.set_ticks_position('none') 
      im_ax.yaxis.set_ticks_position('none') 
      im_ax.tick_params(  axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        left=False,
        right=False,
        labelbottom=False,
        labeltop=False,
        labelleft=False,
        labelright=False )


      im_ax.spines['bottom'].set_color('white')
      im_ax.spines['top'].set_color('white')
      im_ax.spines['left'].set_color('white')
      im_ax.spines['right'].set_color('white')


      im_fig = plt.figure(figsize=(image_size, image_size))
      im_ax = im_fig.add_axes([0.0, 0.0, 1.0, 1.0])
      im_ax.imshow(face_image, cmap='magma')

      im_fig.savefig('./plots/gas_density_image_{}_{}_{}.png'.format(file_index, n_pixels-1, colormap), dpi=n_pixels-1)
      plt.close()

      im_fig = plt.figure(figsize=(image_size, image_size * edge_aspect_ratio))
      im_ax = im_fig.add_axes([0.0, 0.0, 1.0, 1.0])
      im_ax.imshow(edge_image, cmap='magma')

      im_fig.savefig('./plots/gas_density_edge_image_{}_{}_{}.png'.format(file_index, n_pixels-1, colormap), dpi=n_pixels-1)
      plt.close()

#plt.show()

#all_im_fig.subplots_adjust(right=0.875)
#cbar_ax = all_im_fig.add_axes([0.875,0.125,0.125,0.75])
#all_im_fig.colorbar( face_image, cax=cbar_ax, shrink=0.9, label=r'Gas Column Density [${\rm M}_{\odot}$/{\rm pc}^2]')

all_im_fig.savefig('OfficialPaperPlots/F6-GasMap.png', dpi=n_pixels-1)
plt.show()

#      for this_angle in range(90):
#        ls = LightSource(azdeg=315, altdeg=this_angle)
#        this_cmap = matplotlib.colormaps[colormap]
#        shaded_image = ls.shade(image, cmap=this_cmap, blend_mode='hsv')
        
#        im_fig = plt.figure(figsize=(image_size, image_size))
#        im_ax = im_fig.add_axes([0.0, 0.0, 1.0, 1.0])
#        im_ax.imshow(shaded_image) 

#        angle_tag = str(this_angle).zfill(3)
    
#        im_fig.savefig('./plots/gas_density_shaded_image_{}_{}_{}_{}.png'.format(file_index, n_pixels-1, colormap, angle_tag), dpi=n_pixels-1)
#        plt.close()





