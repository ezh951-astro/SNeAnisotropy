import readhaloHDF5
import numpy as np
import readsubfHDF5
import tables
from mpi4py import MPI


base='/n/ghernquist/Illustris/Runs/'
run = 'Illustris-3'
snap = 135
boxsize=75000.0


#=========================================================#
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
n_tasks = size
this_task = rank
#=========================================================#



def write_cf00_galaxy_magnitudes():
    cat = readsubfHDF5.subfind_catalog(base+run+'/output/', snap, keysel=['SubhaloStellarPhotometrics'] )
    cat_mags = np.array(cat.SubhaloStellarPhotometrics)
    n_subs = cat_mags.shape[0]
    PartType=4

    mod_cat_mags = cat_mags[:, :]

    for subnr in np.arange(0,n_subs):
        if subnr % 100 == 0:
            print "loading subhalo "+str(subnr)+" out of "+str(n_subs)

#        print cat_mags[subnr,:]
	if cat_mags[subnr,0] > 1e35:
	    print "no star particle here... move on..."
	else:

            luminosity = np.array(readhaloHDF5.readhalo(base+run+'/output/', 'snap', snap, 'GSPH', PartType, -1, subnr))    
	    # really a magnitude

	    new_luminosity = np.zeros(luminosity.shape)
    	    # band luminosities: U, B, V, K, g, r, i, z
	
	    cf00_dust_factors = np.array([1.461, 1.269, 1.0857, 0.409, 1.129, 0.946, 0.843, 0.763])

#######################################################################
#  convert the formation time to an age in Gyrs 
#######################################################################
	    ft = np.array(readhaloHDF5.readhalo(base+run+'/output/', 'snap', snap, 'GAGE', PartType, -1, subnr))
#	    ft = np.arange(11)/10.0  # use this to test my cosmology age conversion below -- seems to work.

            a0 = 1.0 			### might want to replace with a more generic: fload_box_time(frun,snapnum)
            OmegaLambda = 0.7274
            Omega0      = 0.2726
            factor1 = 2.0 / (3.0 * np.sqrt(OmegaLambda));

            term1 = np.sqrt(OmegaLambda / Omega0) * a0**1.5
            term2 = np.sqrt(1 + OmegaLambda / Omega0 * a0**3.0)
            factor2 = np.log(term1 + term2);

            t0 = factor1 * factor2;

            term1 = np.sqrt(OmegaLambda / Omega0) * ft**1.5
            term2 = np.sqrt(1 + OmegaLambda / Omega0 * ft** 3.0)
            factor2 = np.log(term1 + term2);

            t1 = factor1 * factor2;

            result = t0 - t1;

            hubble            = 3.2407789e-18
            hubbleparam       = 0.7
            sec_per_megayear  = 3.15576e13

            ft = result / (hubble * hubbleparam)#;  #/* now in seconds */
            ft /= sec_per_megayear * 1000#;     /* now in gigayears */
	    ft *= 1e9				# now in years
#######################################################################
#  apply the correction 
#######################################################################
	    young_index = ft < 1e7
	    old_index   = ft > 1e7
	    
	    new_luminosity[young_index, :] = luminosity[young_index, :] + cf00_dust_factors
	    new_luminosity[old_index  , :] = luminosity[old_index  , :] + cf00_dust_factors/3.0

	    luminosity     = 10.0**(-0.4 * luminosity )
	    new_luminosity = 10.0**(-0.4 * new_luminosity)

	    mod_cat_mags[subnr, :] = -2.5 * np.log10( np.sum(new_luminosity , axis=0))
#######################################################################
#  save results 
#######################################################################

    f = tables.openFile("cf00_galaxy_photometrics_"+run+"_"+str(snap)+".hdf5", mode="w")
    f.createArray(f.root, 'CF00_photometrics', mod_cat_mags)
    f.close()
#######################################################################


def find_max_r_of_fuzz():
    cat = readsubfHDF5.subfind_catalog(base+run+'/output/', snap, keysel=['GroupPos', 'GroupMass', "Group_R_Crit200" ] )

    GroupMass    	= np.array(cat.GroupMass)
    n_groups = GroupMass.shape[0]

    PartType = 4

    r_max_array  = np.zeros((n_groups))
    x_min_array  = np.zeros((n_groups))
    x_max_array  = np.zeros((n_groups))
    y_min_array  = np.zeros((n_groups))
    y_max_array  = np.zeros((n_groups))
    z_min_array  = np.zeros((n_groups))
    z_max_array  = np.zeros((n_groups))

    global_r_max_array  = np.zeros((n_groups))
    global_x_min_array  = np.zeros((n_groups))
    global_x_max_array  = np.zeros((n_groups))
    global_y_min_array  = np.zeros((n_groups))
    global_y_max_array  = np.zeros((n_groups))
    global_z_min_array  = np.zeros((n_groups))
    global_z_max_array  = np.zeros((n_groups))

    for groupnr in np.arange(0,n_groups):
      if ((groupnr % n_tasks) == this_task):
	print "processing group "+str(groupnr)+" on task "+str(this_task)


	if groupnr % 100 == 0:
	    print "loading group "+str(groupnr)+" out of "+str(n_groups)
	data = np.array([ [0, 0, 0] ])
	PartTypeArray = [0, 1, 4]
	for PartType in PartTypeArray:
	    tmp = np.array(readhaloHDF5.readhalo(base+run+'/output/', 'snap', snap, 'POS ', PartType, groupnr, -1, load_fuzz=1))
	    if not tmp is None and len(tmp.shape) > 0:
	        data = np.append(data, tmp, axis=0)

	if data is None:
	    print "No data in this group!!!"
        else:
          x = data[:,0] - cat.GroupPos[groupnr,0]
          y = data[:,1] - cat.GroupPos[groupnr,1]
          z = data[:,2] - cat.GroupPos[groupnr,2]

	  x[ x > boxsize/2.0 ] -= boxsize/2.0
          y[ y > boxsize/2.0 ] -= boxsize/2.0
          z[ z > boxsize/2.0 ] -= boxsize/2.0

          x[ x < - boxsize/2.0 ] += boxsize/2.0
          y[ y < - boxsize/2.0 ] += boxsize/2.0
          z[ z < - boxsize/2.0 ] += boxsize/2.0

	  r = np.sqrt( x * x + y * y + z * z)
        
	  xmin = np.min(x) + cat.GroupPos[groupnr,0]
          xmax = np.max(x) + cat.GroupPos[groupnr,0]
          ymin = np.min(y) + cat.GroupPos[groupnr,1]
          ymax = np.max(y) + cat.GroupPos[groupnr,1]
          zmin = np.min(z) + cat.GroupPos[groupnr,2]
          zmax = np.max(z) + cat.GroupPos[groupnr,2]

	  dx = 2.0 * np.max( [x, np.abs(x) ])
          dy = 2.0 * np.max( [y, np.abs(y) ])
          dz = 2.0 * np.max( [z, np.abs(z) ])

          max_r = r.max()

	  r_max_array[groupnr] = max_r
	  x_min_array[groupnr] = xmin
          x_max_array[groupnr] = xmax
          y_min_array[groupnr] = ymin
          y_max_array[groupnr] = ymax
          z_min_array[groupnr] = zmin
          z_max_array[groupnr] = zmax


    comm.Barrier()
    comm.Allreduce(r_max_array, global_r_max_array, op=MPI.SUM)
    comm.Allreduce(x_min_array, global_r_max_array, op=MPI.SUM)
    comm.Allreduce(x_max_array, global_r_max_array, op=MPI.SUM)
    comm.Allreduce(y_min_array, global_r_max_array, op=MPI.SUM)
    comm.Allreduce(y_max_array, global_r_max_array, op=MPI.SUM)
    comm.Allreduce(z_min_array, global_r_max_array, op=MPI.SUM)
    comm.Allreduce(z_max_array, global_r_max_array, op=MPI.SUM)



    f = tables.openFile("maximum_group_fuzz_radii_"+run+"_"+str(snap)+".hdf5", mode="w")
    f.createArray(f.root, 'MaximumR', global_r_max_array)
    f.close()

    f = tables.openFile("bounding_box_for_fuzz_"+run+"_"+str(snap)+".hdf5", mode="w")
    f.createArray(f.root, 'x_min', global_x_min_array)
    f.createArray(f.root, 'x_max', global_x_max_array)
    f.createArray(f.root, 'y_min', global_y_min_array)
    f.createArray(f.root, 'y_max', global_y_max_array)
    f.createArray(f.root, 'z_min', global_z_min_array)
    f.createArray(f.root, 'z_max', global_z_max_array)
    f.close()

    


#n_groups = GroupMass.shape[0]

#for groupnr in np.arange(n_groups):
    










