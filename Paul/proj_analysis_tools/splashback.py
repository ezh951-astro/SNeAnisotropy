# test
import numpy as np
from scipy.signal import savgol_filter
import simread.readsubfHDF5 as rsubf
import simread.cutoutsHDF5 as cut
import util.calc_hsml as calc_hsml

class splashback_analysis:
    def __init__(self, datadir, snapnr, min_radius = 0.1, max_radius = 4, n_radial_bins = 500, n_samples_per_shell = 250, n_ang_bins=100, n_z_bins=100, little_h = 0.6774, **kwargs):
        self.datadir = datadir
        self.snapnr  = snapnr 
        self.little_h = little_h

        self.min_radius = min_radius
        self.max_radius = max_radius

        self.grid_radius = 10.0**np.linspace( np.log10(min_radius), np.log10(max_radius), n_radial_bins, endpoint=True)
    
        self.n_radial_bins = n_radial_bins
        self.n_samples_per_shell = n_samples_per_shell
        self.n_ang_bins = n_ang_bins
        self.n_z_bins   = n_z_bins

        self.set_shell_sampling( ) 
        self.set_cyl_sampling( )
        self.mass_density = np.zeros_like( self.grid_radius )            
        self.density_array = np.zeros_like( self.grid_radius ) 
        self.cat = rsubf.subfind_catalog( self.datadir, self.snapnr , keysel=['GroupPos', 'SubhaloPos', 'GroupFirstSub', 'Group_R_Mean200', 'Group_M_Crit200'] )

        self.sbr_dm_2d_array = np.zeros( len( self.cat.Group_R_Mean200 ) ) 
        self.sbr_dm_3d_array = np.zeros( len( self.cat.Group_R_Mean200 ) ) 

        self.sbr_stellar_2d_array = np.zeros( len( self.cat.Group_R_Mean200 ) ) 
        self.sbr_stellar_3d_array = np.zeros( len( self.cat.Group_R_Mean200 ) ) 



    def set_shell_sampling(self):
        phi = 3.14159 * (3. - np.sqrt(5.))
        points = np.arange(self.n_samples_per_shell)
        self.y_shell = 1 - (1.0*points / (self.n_samples_per_shell-1.0))*2.0
        shell_radius=np.sqrt( 1 - self.y_shell * self.y_shell )
        theta = phi*points
        self.x_shell = np.cos(theta) * shell_radius
        self.z_shell = np.sin(theta) * shell_radius


    def set_cyl_sampling(self):
        theta = np.linspace( 0, 2*3.14159, self.n_ang_bins, endpoint=False )

        self.x_cyl = np.cos(theta)
        self.y_cyl = np.sin(theta)
        self.z_cyl = np.linspace( - self.max_radius  ,  self.max_radius  , self.n_z_bins )
        #note: x/y pairins correspond to points on a circle.  There should be a different circle for each z-point to create the cylinder.


    def assign_fofnr( self, fofnr ):
        self.fofnr       = fofnr
        self.cutout_file = self.datadir+'/cutouts/snapdir_'+str(self.snapnr).zfill(3)+'/fof_box_'+str(fofnr)+'.hdf5'
        self.r200        = self.cat.Group_R_Mean200[fofnr] 

    def load_cutout_data( self ):
        all_data = cut.read_box_cutout_file( self.cutout_file )
        self.pos_dm = all_data["PartType1/Coordinates"]
        self.mass_dm = all_data["PartType1/Masses"] * 1e10 / self.little_h
        self.pos_s  = all_data["PartType4/Coordinates"] 
        self.mass_s = all_data["PartType4/Masses"] * 1e10 / self.little_h
        self.dist_dm = np.sqrt(np.sum( self.pos_dm * self.pos_dm , axis=1 ))
        self.dist_s = np.sqrt(np.sum(  self.pos_s  * self.pos_s , axis=1 ))


    def set_density_array( self, target, mode='3D', target_ngb=32):

        # set the locations for 3D density reconstruction. 
        x_shell_list = np.zeros( self.n_radial_bins * self.n_samples_per_shell )
        y_shell_list = np.zeros( self.n_radial_bins * self.n_samples_per_shell )
        z_shell_list = np.zeros( self.n_radial_bins * self.n_samples_per_shell )
        for ijk in range( len(self.grid_radius) ):                                                          
            start_index = (0 + ijk) * self.n_samples_per_shell 
            end_index   = (1 + ijk) * self.n_samples_per_shell 

            x_shell_list[start_index:end_index] = self.x_shell*self.grid_radius[ijk]*self.r200 
            y_shell_list[start_index:end_index] = self.y_shell*self.grid_radius[ijk]*self.r200     
            z_shell_list[start_index:end_index] = self.z_shell*self.grid_radius[ijk]*self.r200 

        # set the locations for 2D density reconstruction.
        x_cyl_list = np.zeros( self.n_radial_bins * self.n_ang_bins * self.n_z_bins )
        y_cyl_list = np.zeros( self.n_radial_bins * self.n_ang_bins * self.n_z_bins )
        z_cyl_list = np.zeros( self.n_radial_bins * self.n_ang_bins * self.n_z_bins )

        for r_ijk in range( len(self.grid_radius)):
            start_r_index = (r_ijk + 0) * ( self.n_ang_bins * self.n_z_bins )  # this is ijk (the radial bin) times the number of points in each radial bin (n_ang_bins x n_z_bins) 
            end_r_index   = (r_ijk + 1) * ( self.n_ang_bins * self.n_z_bins )
            for phi_ijk in range( self.n_z_bins ):
                start_phi_index = start_r_index + (phi_ijk+0) * self.n_z_bins
                end_phi_index   = start_r_index + (phi_ijk+1) * self.n_z_bins
 
                x_cyl_list[start_phi_index:end_phi_index] = self.x_cyl[phi_ijk]*self.grid_radius[r_ijk]*self.r200   # all the same x values
                y_cyl_list[start_phi_index:end_phi_index] = self.y_cyl[phi_ijk]*self.grid_radius[r_ijk]*self.r200   # all the same y values
                z_cyl_list[start_phi_index:end_phi_index] = self.z_cyl * self.r200 

        print("calculating density on sphere points...")
        if target=='stars':
            density_list_3d = calc_hsml.get_gas_density_around_stars( self.pos_s[:,0],  self.pos_s[:,1],  self.pos_s[:,2],  self.mass_s,  x_shell_list, y_shell_list, z_shell_list, DesNgb=target_ngb)
            density_list_2d = calc_hsml.get_gas_density_around_stars( self.pos_s[:,0],  self.pos_s[:,1],  self.pos_s[:,2],  self.mass_s,  x_cyl_list,   y_cyl_list,   z_cyl_list,   DesNgb=target_ngb)
        elif target=='DM':
            density_list_3d = calc_hsml.get_gas_density_around_stars( self.pos_dm[:,0], self.pos_dm[:,1], self.pos_dm[:,2], self.mass_dm, x_shell_list, y_shell_list, z_shell_list, DesNgb=target_ngb)
            density_list_2d = calc_hsml.get_gas_density_around_stars( self.pos_dm[:,0], self.pos_dm[:,1], self.pos_dm[:,2], self.mass_dm, x_cyl_list,   y_cyl_list,   z_cyl_list,   DesNgb=target_ngb)
    
        density_array_3d = np.zeros( self.n_radial_bins )
        density_array_2d = np.zeros( self.n_radial_bins )

        #collapse the 2d densities into surface densities
        for r_ijk in range(self.n_radial_bins):
            this_density_array = np.zeros( self.n_ang_bins )
            for phi_ijk in range(self.n_ang_bins):
                start_index = r_ijk*(self.n_ang_bins*self.n_z_bins) + (phi_ijk+0)*self.n_z_bins
                end_index   = r_ijk*(self.n_ang_bins*self.n_z_bins) + (phi_ijk+1)*self.n_z_bins
                this_density_array[phi_ijk] = np.sum( density_list_2d[ start_index : end_index ] ) * (np.max(self.z_cyl) - np.min(self.z_cyl) ) * self.r200 / self.n_z_bins # for each phi, need so integrate alone LOS
            density_array_2d[r_ijk] = np.median( this_density_array )

        for r_ijk in range(self.n_radial_bins):
            this_density_list = density_list_3d[(r_ijk * self.n_samples_per_shell):((r_ijk+1)*self.n_samples_per_shell) ]
            density_array_3d[r_ijk] = np.median( this_density_list )

    
        if target=='stars':
          self.stellar_density_array = density_array_3d
          self.stellar_surf_density_array = density_array_2d
        elif target=='DM':
          self.dm_density_array      =    density_array_3d
          self.dm_surf_density_array      = density_array_2d
    
    
    
    def set_gradient_array( self, target, smooth_window_density=5, smooth_order_density=1, smooth_window_grad=31, smooth_order_grad=1 ):
        this_log_radius = np.log10( self.grid_radius )
        if target=='stars':
            log_density_3d     = np.log10( self.stellar_density_array )
            log_density_2d     = np.log10( self.stellar_surf_density_array ) 
        elif target=='DM':
            log_density_3d     = np.log10( self.dm_density_array ) 
            log_density_2d     = np.log10( self.dm_surf_density_array )
    
        log_density_2d = savgol_filter(log_density_2d,smooth_window_density,smooth_order_density)    # 31,4) ## 31, 4)
        log_density_3d = savgol_filter(log_density_3d,smooth_window_density,smooth_order_density)    # 31,4) ## 31, 4)
    
        gradient_2d = savgol_filter( np.gradient( log_density_2d, this_log_radius, edge_order=2) , smooth_window_grad, smooth_order_grad) # this is dependent on the number of radial bins you are using...
        gradient_3d = savgol_filter( np.gradient( log_density_3d, this_log_radius, edge_order=2) , smooth_window_grad, smooth_order_grad) # this is dependent on the number of radial bins you are using...


        if target=='stars':
            self.smoothed_stellar_density_array       = 10.0**log_density_3d
            self.smoothed_stellar_surf_density_array  = 10.0**log_density_2d
            self.smoothed_stellar_gradient_array      = gradient_3d
            self.smoothed_stellar_surf_gradient_array = gradient_2d
        elif target=='DM':
            self.smoothed_dm_density_array            = 10.0**log_density_3d
            self.smoothed_dm_surf_density_array       = 10.0**log_density_2d
            self.smoothed_dm_gradient_array           = gradient_3d
            self.smoothed_dm_surf_gradient_array      = gradient_2d
 
    def find_sbr( self, target, mode='3D' ):
        this_radius = self.grid_radius 
        if target=='stars':
            if mode=='3D':
                this_gradient = self.smoothed_stellar_gradient_array
            elif mode=='2D':
                this_gradient = self.smoothed_stellar_surf_gradient_array
        elif target=='DM':
            if mode=='3D':
                this_gradient = self.smoothed_dm_gradient_array
            elif mode=='2D':
                this_gradient = self.smoothed_dm_surf_gradient_array

        index_array = np.array( range( len( this_gradient) )  ) 
        print( np.min(this_gradient) ) 
        print( this_gradient == np.min(this_gradient) ) 
        print( index_array[ this_gradient == np.min(this_gradient) ] ) 
        min_index = index_array[ this_gradient == np.min(this_gradient) ]

        print(min_index)
      
        while(min_index==0 or (min_index == len(this_gradient)-1) ): 
            if(min_index==0):
                this_gradient=this_gradient[1:]
                this_radius  = this_radius[1:]
            else:
                this_gradient=this_gradient[:-2]
                this_radius  = this_radius[:-2]
            min_index = index_array[ this_gradient == np.min(this_gradient) ]


        if target=='stars':
            if mode=='3D':
                self.sbr_stellar_3d_array[ self.fofnr ] = this_radius[min_index] 
            elif mode=='2D':
                self.sbr_stellar_2d_array[ self.fofnr ] = this_radius[min_index] 
        elif target=='DM':
            if mode=='3D':
                self.sbr_dm_3d_array[ self.fofnr ] = this_radius[min_index] 
            elif mode=='2D':
                self.sbr_dm_2d_array[ self.fofnr ] = this_radius[min_index] 

        return this_radius[min_index]


    def return_radius( self ):
        return self.grid_radius

    def return_density( self, target, mode='3D' ):
        if mode=='3D':
            if target=='stars':
                return self.smoothed_stellar_density_array
            elif target=='DM':
                return self.smoothed_dm_density_array
        elif mode=='2D':
            if target=='stars':
                return self.smoothed_stellar_surf_density_array
            elif target=='DM':
                return self.smoothed_dm_surf_density_array


    def return_gradient( self, target, mode='3D' ):
        if mode=='3D':
            if target=='stars':
                return self.smoothed_stellar_gradient_array
            elif target=='DM':
                return self.smoothed_dm_gradient_array + 7.5
        elif mode=='2D':
            if target=='stars':
                return self.smoothed_stellar_surf_gradient_array
            elif target=='DM':
                return self.smoothed_dm_surf_gradient_array + 7.5
 



def find_matching_fofs( fofnr, splashback_1, splashback_2, splashback_3 ):
    target_pos = splashback_3.cat.GroupPos[fofnr,:]

    offsets_pos_1 = splashback_1.cat.GroupPos[:,:] - target_pos
    offsets_pos_2 = splashback_2.cat.GroupPos[:,:] - target_pos
    offsets_pos_3 = splashback_3.cat.GroupPos[:,:] - target_pos

    offsets_1 = np.sqrt( offsets_pos_1[:,0]*offsets_pos_1[:,0] + offsets_pos_1[:,1]*offsets_pos_1[:,1] + offsets_pos_1[:,2]*offsets_pos_1[:,2] ) / splashback_1.cat.Group_R_Mean200[:]
    offsets_2 = np.sqrt( offsets_pos_2[:,0]*offsets_pos_2[:,0] + offsets_pos_2[:,1]*offsets_pos_2[:,1] + offsets_pos_2[:,2]*offsets_pos_2[:,2] ) / splashback_2.cat.Group_R_Mean200[:]
    offsets_3 = np.sqrt( offsets_pos_3[:,0]*offsets_pos_3[:,0] + offsets_pos_3[:,1]*offsets_pos_3[:,1] + offsets_pos_3[:,2]*offsets_pos_3[:,2] ) / splashback_3.cat.Group_R_Mean200[:]

    fofnr_1 = np.argmin( offsets_1[:20000])
    fofnr_2 = np.argmin( offsets_2[:20000])
    fofnr_3 = np.argmin( offsets_3[:20000])

    print("\n\n\n")
    print("Closest FoF in TNG300-1 is Fofnr {} which is {} Mpc away which is {:.2} of the virial radius".format( fofnr_3, offsets_3[fofnr_3], offsets_3[fofnr_3] / splashback_3.cat.Group_R_Mean200[fofnr_3] ) )
    print("Closest FoF in TNG300-2 is Fofnr {} which is {} Mpc away which is {:.2} of the virial radius".format( fofnr_2, offsets_2[fofnr_2], offsets_2[fofnr_2] / splashback_2.cat.Group_R_Mean200[fofnr_2] ) )
    print("Closest FoF in TNG300-3 is Fofnr {} which is {} Mpc away which is {:.2} of the virial radius".format( fofnr_1, offsets_1[fofnr_1], offsets_1[fofnr_1] / splashback_1.cat.Group_R_Mean200[fofnr_1] ) )

    return fofnr_1, fofnr_2, fofnr_3

