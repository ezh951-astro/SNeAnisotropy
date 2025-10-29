#!/usr/bin/env python
""" routines for reading and writing snapshot cutout data -- likely best suited for cosmo sims.

# Python HDF5 snapshot reader/writer for AREPO/GADGET
#
# import snapHDF5 as snap
# header = snap.snapshot_header("snap_063")
# mass = snap.read_block("snap_063", "MASS", parttype=5)
# pos = snap.read_block("snap_063", "POS ", parttype=1)
# printpos
#
#
"""

__author__ = "Paul Torrey"
__copyright__ = "Copyright 2022, The Authors"
__credits__ = ["Paul Torrey"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Paul Torrey"
__email__ = "paul.torrey@ufl.edu"

import util.naming as naming
import h5py
import simread.readsnapHDF5 as rs
import numpy as np


def make_box_cutout_file( datadir, snapnr, cutout_file, cutout_center, cutout_size, boxsize=None):

    filenames = naming.get_snap_filenames( datadir, 'snapdir', snapnr)
    filenames.sort()

    with h5py.File(cutout_file, "w") as f:
      for parttype in [0,1,4,5]:
        pos = None
        mass = None
        vel  = None
        for filenr,name in enumerate(filenames):
          print(name )
          this_pos = rs.read_block_single_file(name, 'Coordinates', 3, parttype=parttype)[0]

          for ijk in range(3):
            this_pos[:,ijk] -= cutout_center[ijk]
            if boxsize is not None:
              this_pos[this_pos[:,ijk] >  boxsize/2.0 ] -= boxsize
              this_pos[this_pos[:,ijk] < -boxsize/2.0 ] += boxsize

          in_box = ( this_pos[:,0] > -cutout_size ) & ( this_pos[:,0] < cutout_size ) & \
                   ( this_pos[:,1] > -cutout_size ) & ( this_pos[:,1] < cutout_size ) & \
                   ( this_pos[:,2] > -cutout_size ) & ( this_pos[:,2] < cutout_size )

          if np.sum(in_box) > 0: # we have at least one particle in the box on this file
            this_mass =  rs.read_block_single_file(name, 'Masses', 1, parttype=parttype)[0]
            this_vel  =  rs.read_block_single_file(name, 'Velocities', 3, parttype=parttype)[0]
 
            if parttype==4:
              this_gfm_lum = rs.read_block_single_file(name, 'GFM_StellarPhotometrics', 8, parttype=parttype)[0]
              this_gfm_formation_time = rs.read_block_single_file(name, 'GFM_StellarFormationTime', 1, parttype=parttype)[0]
              this_gfm_metallicity = rs.read_block_single_file(name, 'GFM_Metallicity', 1, parttype=parttype)[0]
              this_gfm_initial_mass = rs.read_block_single_file(name, 'GFM_InitialMass', 1, parttype=parttype)[0]

            if pos is None:
                pos = np.array(  this_pos[in_box,:] )
                mass = np.array( this_mass[in_box] )
                vel  = np.array( this_vel[in_box,:] )
                if parttype==4:
                    lum = np.array( this_gfm_lum[in_box,:] )
                    ft  = np.array( this_gfm_formation_time[in_box] )
                    met = np.array( this_gfm_metallicity[in_box] )
                    imass = np.array( this_gfm_initial_mass[in_box] )
            else:
                pos = np.append( pos,  this_pos[in_box,:] ,axis=0 )
                mass= np.append( mass, this_mass[in_box] )
                vel = np.append( vel,  this_vel[in_box,:] ,axis=0 )
                if parttype==4:
                    lum   = np.append( lum,   this_gfm_lum[in_box,:], axis=0 )
                    ft    = np.append( ft,    this_gfm_formation_time[in_box] )
                    met   = np.append( met,   this_gfm_metallicity[in_box] )
                    imass = np.append( imass, this_gfm_initial_mass[in_box] )


        print("All particles for parttype {:d} are loaded".format( parttype ) )
        print("  The min/max x particle positions are {:16.8f}/{:16.8f}".format( np.min(pos[:,0]), np.max(pos[:,0]) ) )
        print("  The min/max y particle positions are {:16.8f}/{:16.8f}".format( np.min(pos[:,1]), np.max(pos[:,1]) ) )
        print("  The min/max z particle positions are {:16.8f}/{:16.8f}".format( np.min(pos[:,2]), np.max(pos[:,2]) ) )

        pt = f.create_group("PartType"+str(parttype))
        pt.create_dataset( "Coordinates" , data = pos )
        pt.create_dataset( "Velocities" , data = vel )
        pt.create_dataset( "Masses" , data = mass )
        if parttype==4:
            pt.create_dataset( "GFM_StellarPhotometrics", data = lum )
            pt.create_dataset( "GFM_StellarFormationTime", data = ft )
            pt.create_dataset( "GFM_Metallicity", data = met )
            pt.create_dataset( "GFM_InitialMass", data = imass )





def read_box_cutout_file( cutout_file ):
  data = {}
  with h5py.File( cutout_file, 'r' ) as f:
      for PartType in f:
          for field in f[PartType]:
              tag = "{}/{}".format( PartType, field ) 
              data[tag] = np.array( f[PartType][field][:] ) 

  return data 

