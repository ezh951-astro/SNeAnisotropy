import numpy as np
import h5py
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
zdirsname = 'zdirs'

### handle the input file directories

zdirslist = open(zdirsname,'r').read().splitlines()

header_indices = []
lineskip_indices = []
for i,linestr in enumerate(zdirslist):
	if linestr == '':
		lineskip_indices.append(i)
	elif linestr[0] == '%':
		header_indices.append(i)

Nruns = lineskip_indices[0] - header_indices[0] - 1

runlabels      = zdirslist[ header_indices[0]+1 : header_indices[0]+Nruns+1 ]
evodata_dirs   = zdirslist[ header_indices[1]+1 : header_indices[1]+Nruns+1 ]
finalsnap_dirs = zdirslist[ header_indices[2]+1 : header_indices[2]+Nruns+1 ]
finalfof_dirs  = zdirslist[ header_indices[3]+1 : header_indices[3]+Nruns+1 ]

if evodata_dirs[0] == 'None':
    DEFAULT_MODE = 'ISO'
    print('Evodata None Detected; runs are ISO')
else:
    DEFAULT_MODE = 'COSMO'
    print('Evodata Loaded; runs are COSMO')

snap_dirs = [ dirstr[:-8] for dirstr in finalsnap_dirs ]
fof_dirs = [ dirstr[:-8] for dirstr in finalfof_dirs ]

final_snap_nums = [ int(dirstr[-8:-5]) for dirstr in finalsnap_dirs ]
print( dirstr[:-3] for dirstr in finalsnap_dirs )

# both ISO and COSMO have snapshots:
finalsnap_list = [ h5py.File(fname) for fname in finalsnap_dirs ]

print(runlabels)
print(finalsnap_list)

# only COSMO has evodata and finalfof.
if DEFAULT_MODE == 'COSMO':
    evodata_list   = [ h5py.File(fname) for fname in evodata_dirs   ]
    finalfof_list  = [ h5py.File(fname) for fname in finalfof_dirs  ]
    print(evodata_list)
    print(finalfof_list)


### define some global constants

G = 4.302e-6 # kpc Msun^-1 (km/s)^2
ageGyr = 13.77 # age of universe in Gyr

### things about plot styles

linestyles = ['solid','dashed','dotted']
linecolors_global = ['red','blue','green']

##### FUNCTIONS

def handle_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--RunNums', help='which run (numbers) do we want to include')
    parser.add_argument('-v', '--version', help='string noting the version, use todays date')
    args = parser.parse_args()

    if args.RunNums != None:
        ActiveRunIDs = [ int(args.RunNums[i]) for i in range(len(args.RunNums)) ]
    else:
        ActiveRunIDs = range(Nruns)

    if args.version != None:
        date = args.version
    else:
        date = 'xxxx'

    return ActiveRunIDs, date

def get_istr3(i):
    # Assumes that i>=0
    if i<10:
        return '00'+str(i)
    if i>=10 and i<100:
        return '0'+str(i)
    if i>=100:
        return str(i)

def create_cmap(cmapname,lo=4.0,hi=8.0,Nbins=20):
    levels = MaxNLocator(nbins=Nbins).tick_values(lo,hi)
    cmap = plt.get_cmap(cmapname)
    norm = BoundaryNorm(levels,ncolors=cmap.N,clip=True)
    return cmap, norm

def order_masses( rList, mList ):
    ind = np.argsort(rList)
    rOrd = rList[ind]
    mOrd = mList[ind]
    cummass = np.cumsum(mOrd)
    return rOrd, cummass

def half_mass_radius( rList, mList ):

    ordered = order_masses(rList,mList)

    distOrd = ordered[0]
    cummass = ordered[1]
    #print(distOrd,cummass)

    TotMass = cummass[-1]
    HalfMass = TotMass / 2.
    indHalf = np.argmax( cummass>HalfMass )
    rHalf = distOrd[indHalf]

    return rHalf

def GasTemp(f,mask0=None):

    # Temperature of Gas Particles of a given Snapshot

    UnitLength_in_cm =        3.085678e21 #  1.0 kpc, but this doesn't matter
    UnitMass_in_g    =        1.989e43 #  1.0e10 solar masses 
    UnitVelocity_in_cm_per_s = 1.e5  #  1 km/sec 
    UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s
    #UnitDensity_in_cgs= UnitMass_in_g/ UnitLength_in_cm**3
    #UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ UnitTime_in_s**2
    UnitEnergy_in_cgs= UnitMass_in_g * UnitLength_in_cm**2 / UnitTime_in_s**2

    BOLTZMANN = 1.3806e-16
    PROTONMASS = 1.6726e-24

    Xh=0.76                        # mass fraction of hydrogen
    gamma= 5.0/3.

    U = np.array(f[u'PartType0'][u'InternalEnergy'])
    Nelec = np.array(f[u'PartType0'][u'ElectronAbundance'])

    if mask0 is None:
      mask0 = (U>0.)

    U = U[mask0]
    Nelec = Nelec[mask0]

    MeanWeight= 4.0/(3*Xh+1+4*Xh*Nelec) * PROTONMASS
    temp = MeanWeight/BOLTZMANN * (gamma-1) * U * UnitEnergy_in_cgs/ UnitMass_in_g

    return temp

def MainSubhaloVelocity(f,s,N=128):

    # CM velocity of the N=128 closest particles to the halo center. Returns it in PHYSICAL units

    redshift = f[u'Header'].attrs[u'Redshift']
    scalefac = 1. / (1.+redshift)
    h_small = f[u'Header'].attrs['HubbleParam']
    SimTime = f[u'Header'].attrs[u'Time']

    center = np.array(s[u'Subhalo'][u'SubhaloPos'][0]) * 1000.*scalefac/h_small

    xyz0 = np.array(f[u'PartType0'][u'Coordinates'])*1000.*scalefac/h_small - center
    xyz1 = np.array(f[u'PartType1'][u'Coordinates'])*1000.*scalefac/h_small - center
    xyz4 = np.array(f[u'PartType4'][u'Coordinates'])*1000.*scalefac/h_small - center

    dist0 = np.sqrt( xyz0[:,0]*xyz0[:,0] + xyz0[:,1]*xyz0[:,1] + xyz0[:,2]*xyz0[:,2] )
    dist1 = np.sqrt( xyz1[:,0]*xyz1[:,0] + xyz1[:,1]*xyz1[:,1] + xyz1[:,2]*xyz1[:,2] )
    dist4 = np.sqrt( xyz4[:,0]*xyz4[:,0] + xyz4[:,1]*xyz4[:,1] + xyz4[:,2]*xyz4[:,2] )

    v0 = np.array(f[u'PartType0'][u'Velocities'])*np.sqrt(scalefac)/h_small #units are now physical km/s
    v1 = np.array(f[u'PartType1'][u'Velocities'])*np.sqrt(scalefac)/h_small
    v4 = np.array(f[u'PartType4'][u'Velocities'])*np.sqrt(scalefac)/h_small
    
    mass0 = np.array(f[u'PartType0'][u'Masses'])*10**10/h_small
    mass1 = np.array(f[u'Header'].attrs[u'MassTable'][1])*10**10/h_small * np.ones(np.shape(xyz1)[0])
    mass4 = np.array(f[u'PartType4'][u'Masses'])*10**10/h_small

    dist = np.concatenate((dist0,dist1,dist4))
    mass = np.concatenate((mass0,mass1,mass4))
    vel  = np.concatenate((v0,v1,v4))
    
    ind = np.argsort(dist)
    
    massNbhd = mass[ind][:N]
    velNbhd  = vel[ind][:N]
    
    Mtot = sum(massNbhd)
    velWeight = [ massNbhd[i]*velNbhd[i] for i in range(N) ]
    #print(np.shape(velNbhd))
    CMvel = sum(velWeight) / Mtot
    
    return CMvel

def ptcl_density_bin(r1,r2,mass,min1=-50.,max1=-50.,min2=-50.,max2=50.,Nbins1=100,Nbins2=100):

    ptcl_density = np.zeros((Nbins1,Nbins2))

    for i_ptcl, (x,y,m) in enumerate(zip(r1,r2,mass)):
        xbin = int( (x-min1) * Nbins1 / (max1-min1) )
        ybin = int( (y-min2) * Nbins2 / (max2-min2) )
        if xbin>=0 and ybin>=0:
            try:
                ptcl_density[xbin,ybin] += m
            except IndexError:
                pass
        else:
            pass
        
    ptcl_density = ptcl_density * Nbins1 * Nbins2 / (max1-min1) / (max2-min2) # divide by area of bin
    
    return np.transpose(ptcl_density)

def morphology(f,mode=DEFAULT_MODE,origin=None,s=None,rad=50.,nbins=100,ptypes=[0,1,4],mask0=None,mask1=None,mask4=None):

    proj0_all = np.zeros((3,nbins,nbins))
    proj1_all = np.zeros((3,nbins,nbins))
    proj4_all = np.zeros((3,nbins,nbins))

    if mode=='ISO':
        redshift = 0.
        scalefac = 1.
        h_small = 1.
        dist_unit_fac = 1.
        origin = np.array(origin)

    if mode=='COSMO':
        redshift = f[u'Header'].attrs[u'Redshift']
        scalefac = 1. / (1.+redshift)
        h_small = f[u'Header'].attrs['HubbleParam']
        dist_unit_fac = 1000.
        origin = np.array(s[u'Subhalo'][u'SubhaloPos'][0]) * dist_unit_fac * scalefac / h_small

    if 0 in ptypes:

        xyz0 = np.array(f[u'PartType0'][u'Coordinates'])*dist_unit_fac*scalefac/h_small - origin
        mass0 = np.array(f[u'PartType0'][u'Masses'])*1.e10/h_small

        if mask0 is None:
            mask0 = (mass0>0.)
        kk0 = (np.abs(xyz0[:,0])<rad)*(np.abs(xyz0[:,1])<rad)*(np.abs(xyz0[:,2])<rad) * mask0
        xyz0 = xyz0[kk0]
        mass0 = mass0[kk0]
        dist0 = np.sqrt( xyz0[:,0]*xyz0[:,0] + xyz0[:,1]*xyz0[:,1] + xyz0[:,2]*xyz0[:,2] )

        proj0_all[0] = ptcl_density_bin( xyz0[:,0], xyz0[:,1], mass0 , min1=-rad,max1=rad,min2=-rad,max2=rad,Nbins1=nbins,Nbins2=nbins ) #xy
        proj0_all[1] = ptcl_density_bin( xyz0[:,0], xyz0[:,2], mass0 , min1=-rad,max1=rad,min2=-rad,max2=rad,Nbins1=nbins,Nbins2=nbins ) #xz
        proj0_all[2] = ptcl_density_bin( xyz0[:,1], xyz0[:,2], mass0 , min1=-rad,max1=rad,min2=-rad,max2=rad,Nbins1=nbins,Nbins2=nbins ) #yz


    if 1 in ptypes:

        xyz1 = np.array(f[u'PartType1'][u'Coordinates'])*dist_unit_fac*scalefac/h_small - origin
        mass1 = np.array(f[u'Header'].attrs[u'MassTable'][1])*1.e10/h_small * np.ones(np.shape(xyz1)[0])

        if mask1 is None:
            mask1 = (mass1>0.)
        kk1 = (np.abs(xyz1[:,0])<rad)*(np.abs(xyz1[:,1])<rad)*(np.abs(xyz1[:,2])<rad) * mask1
        xyz1 = xyz1[kk1]
        mass1 = mass1[kk1]
        dist1 = np.sqrt( xyz1[:,0]*xyz1[:,0] + xyz1[:,1]*xyz1[:,1] + xyz1[:,2]*xyz1[:,2] )

        proj1_all[0] = ptcl_density_bin( xyz1[:,0], xyz1[:,1], mass1 , min1=-rad,max1=rad,min2=-rad,max2=rad,Nbins1=nbins,Nbins2=nbins )
        proj1_all[1] = ptcl_density_bin( xyz1[:,0], xyz1[:,2], mass1 , min1=-rad,max1=rad,min2=-rad,max2=rad,Nbins1=nbins,Nbins2=nbins )
        proj1_all[2] = ptcl_density_bin( xyz1[:,1], xyz1[:,2], mass1 , min1=-rad,max1=rad,min2=-rad,max2=rad,Nbins1=nbins,Nbins2=nbins )


    if 4 in ptypes:

        xyz4 = np.array(f[u'PartType4'][u'Coordinates'])*dist_unit_fac*scalefac/h_small - origin
        mass4 = np.array(f[u'PartType4'][u'Masses'])*1.e10/h_small

        if mask4 is None:
            mask4 = (mass4>0.)
        kk4 = (np.abs(xyz4[:,0])<rad)*(np.abs(xyz4[:,1])<rad)*(np.abs(xyz4[:,2])<rad) * mask4
        xyz4 = xyz4[kk4]
        mass4 = mass4[kk4]
        dist4 = np.sqrt( xyz4[:,0]*xyz4[:,0] + xyz4[:,1]*xyz4[:,1] + xyz4[:,2]*xyz4[:,2] )

        proj4_all[0] = ptcl_density_bin( xyz4[:,0], xyz4[:,1], mass4 , min1=-rad,max1=rad,min2=-rad,max2=rad,Nbins1=nbins,Nbins2=nbins )
        proj4_all[1] = ptcl_density_bin( xyz4[:,0], xyz4[:,2], mass4 , min1=-rad,max1=rad,min2=-rad,max2=rad,Nbins1=nbins,Nbins2=nbins )
        proj4_all[2] = ptcl_density_bin( xyz4[:,1], xyz4[:,2], mass4 , min1=-rad,max1=rad,min2=-rad,max2=rad,Nbins1=nbins,Nbins2=nbins )

    return proj0_all, proj1_all, 'bigwhale', 'bigwhale', proj4_all

def velocity_phase_data(f,mode=DEFAULT_MODE,origin=None,s=None,ptypes=[0,1,4],mask0=None,mask1=None,mask4=None):

    if mode=='ISO':
        redshift = 0.
        scalefac = 1.
        h_small = 1.
        dist_unit_fac = 1.
        origin = np.array(origin)
        center_velocity = np.array([0,0,0])

    if mode=='COSMO':
        redshift = f[u'Header'].attrs[u'Redshift']
        scalefac = 1. / (1.+redshift)
        h_small = f[u'Header'].attrs['HubbleParam']
        dist_unit_fac = 1000.
        origin = np.array(s[u'Subhalo'][u'SubhaloPos'][0]) * dist_unit_fac * scalefac / h_small
        center_velocity = MainSubhaloVelocity(f,s)    

    if 0 in ptypes:

        xyz0 = np.array(f[u'PartType0'][u'Coordinates'])*dist_unit_fac*scalefac/h_small - origin
        v0 = np.array(f[u'PartType0'][u'Velocities'])*np.sqrt(scalefac)/h_small - center_velocity
        mass0 = np.array(f[u'PartType0'][u'Masses'])*10**10/h_small
        dist0 = np.sqrt( xyz0[:,0]*xyz0[:,0] + xyz0[:,1]*xyz0[:,1] + xyz0[:,2]*xyz0[:,2] )
        speed0 = np.sqrt( v0[:,0]*v0[:,0] + v0[:,1]*v0[:,1] + v0[:,2]*v0[:,2] )
        radvel0 = ( xyz0[:,0]*v0[:,0] + xyz0[:,1]*v0[:,1] + xyz0[:,2]*v0[:,2] ) / (dist0+0.001)

        if mask0 is None:
            mask0 = (mass0>0.)

        mass0 = mass0[mask0]
        dist0 = dist0[mask0]
        speed0 = speed0[mask0]
        radvel0 = radvel0[mask0]

        gasPhaseData = np.array([ mass0, dist0, speed0, radvel0 ])

    else:
        gasPhaseData = np.zeros((4,2))

    if 1 in ptypes:

        xyz1 = np.array(f[u'PartType1'][u'Coordinates'])*dist_unit_fac*scalefac/h_small - origin
        v1 = np.array(f[u'PartType1'][u'Velocities'])*np.sqrt(scalefac)/h_small - center_velocity
        mass1 = np.array(f[u'Header'].attrs[u'MassTable'][1])*10**10/h_small * np.ones(np.shape(xyz1)[0])
        dist1 = np.sqrt( xyz1[:,0]*xyz1[:,0] + xyz1[:,1]*xyz1[:,1] + xyz1[:,2]*xyz1[:,2] )
        speed1 = np.sqrt( v1[:,0]*v1[:,0] + v1[:,1]*v1[:,1] + v1[:,2]*v1[:,2] )
        radvel1 = ( xyz1[:,0]*v1[:,0] + xyz1[:,1]*v1[:,1] + xyz1[:,2]*v1[:,2] ) / (dist1+0.001)

        if mask1 is None:
            mask1 = (mass1>0.)    
    
        mass1 = mass1[mask1]
        dist1 = dist1[mask1]
        speed1 = speed1[mask1]
        radvel1 = radvel1[mask1]

        dmPhaseData = np.array([ mass1, dist1, speed1, radvel1 ])

    else:
        dmPhaseData = np.zeros((4,2))

    if 4 in ptypes:

        xyz4 = np.array(f[u'PartType4'][u'Coordinates'])*dist_unit_fac*scalefac/h_small - origin
        v4 = np.array(f[u'PartType4'][u'Velocities'])*np.sqrt(scalefac)/h_small - center_velocity
        mass4 = np.array(f[u'PartType4'][u'Masses'])*10**10/h_small
        dist4 = np.sqrt( xyz4[:,0]*xyz4[:,0] + xyz4[:,1]*xyz4[:,1] + xyz4[:,2]*xyz4[:,2] )
        speed4 = np.sqrt( v4[:,0]*v4[:,0] + v4[:,1]*v4[:,1] + v4[:,2]*v4[:,2] )
        radvel4 = ( xyz4[:,0]*v4[:,0] + xyz4[:,1]*v4[:,1] + xyz4[:,2]*v4[:,2] ) / (dist4+0.001)

        if mask4 is None:
            mask4 = (mass4>0.)
        
        mass4 = mass4[mask4]
        dist4 = dist4[mask4]
        speed4 = speed4[mask4]
        radvel4 = radvel4[mask4]

        starPhaseData = np.array([ mass4, dist4, speed4, radvel4 ])

    else:
        starPhaseData = np.zeros((4,2))

    return gasPhaseData, dmPhaseData, 'bigwhale', 'bigwhale', starPhaseData

def kkRad(f,ptype,mode=DEFAULT_MODE,origin=None,s=None,rad=10.):

    if mode=='ISO':
        redshift = 0.
        scalefac = 1.
        h_small = 1.
        dist_unit_fac = 1.
        origin = np.array(origin)

    if mode=='COSMO':
        redshift = f[u'Header'].attrs[u'Redshift']
        scalefac = 1. / (1.+redshift)
        h_small = f[u'Header'].attrs['HubbleParam']
        dist_unit_fac = 1000.
        origin = np.array(s[u'Subhalo'][u'SubhaloPos'][0]) * dist_unit_fac * scalefac / h_small  

    if ptype==0:
        xyz0 = np.array(f[u'PartType0'][u'Coordinates'])*dist_unit_fac*scalefac/h_small - origin
        dist0 = np.sqrt( xyz0[:,0]*xyz0[:,0] + xyz0[:,1]*xyz0[:,1] + xyz0[:,2]*xyz0[:,2] )
        kk = dist0<rad

    if ptype==1:
        xyz1 = np.array(f[u'PartType1'][u'Coordinates'])*dist_unit_fac*scalefac/h_small - origin
        dist1 = np.sqrt( xyz1[:,0]*xyz1[:,0] + xyz1[:,1]*xyz1[:,1] + xyz1[:,2]*xyz1[:,2] )
        kk = dist1<rad

    if ptype==4:
        xyz4 = np.array(f[u'PartType4'][u'Coordinates'])*dist_unit_fac*scalefac/h_small - origin
        dist4 = np.sqrt( xyz4[:,0]*xyz4[:,0] + xyz4[:,1]*xyz4[:,1] + xyz4[:,2]*xyz4[:,2] )
        kk = dist4<rad
    
    return kk

def kkNewStars(f,dt=50.,mode=DEFAULT_MODE):
    
    SimTime = f[u'Header'].attrs[u'Time']
    StarBirth = np.array(f[u'PartType4'][u'GFM_StellarFormationTime'])

    if mode=='ISO':
        kk4 = StarBirth > SimTime-(dt/1000.)
    if mode=='COSMO':
        kk4 = StarBirth > SimTime-(dt/ageGyr/1000.)
    
    return kk4

def mass_distribution_spherical(f,radArr,mode=DEFAULT_MODE,origin=None,s=None,ptypes=[0,1,4]):

    if mode=='ISO':
        redshift = 0.
        scalefac = 1.
        h_small = 1.
        dist_unit_fac = 1.
        origin = np.array(origin)

    if mode=='COSMO':
        redshift = f[u'Header'].attrs[u'Redshift']
        scalefac = 1. / (1.+redshift)
        h_small = f[u'Header'].attrs['HubbleParam']
        dist_unit_fac = 1000.
        origin = np.array(s[u'Subhalo'][u'SubhaloPos'][0]) * dist_unit_fac * scalefac / h_small

    N = len(radArr)
    ans = np.zeros((7,4,N)) #[ptype][dist,mass,dens,vel][rad]

    if 0 in ptypes:

        xyz0 = np.array(f[u'PartType0'][u'Coordinates']) * dist_unit_fac*scalefac/h_small - origin
        m0 = np.array(f[u'PartType0'][u'Masses']) * 10**10/h_small
        dist0 = np.sqrt( xyz0[:,0]*xyz0[:,0] + xyz0[:,1]*xyz0[:,1] + xyz0[:,2]*xyz0[:,2] )

        m0 = np.concatenate(( m0, 1.e-3*np.ones_like(radArr), 1.e-3*np.ones_like(radArr) ))
        dist0 = np.concatenate(( dist0, radArr+1.e-6, radArr-1.e6 )) # add teeny tiny particles just inside and outside each radius we want to measure at

    else:
        m0 = np.zeros_like(radArr)
        dist0 = np.zeros_like(radArr)


    if 1 in ptypes:

        xyz1 = np.array(f[u'PartType1'][u'Coordinates']) * dist_unit_fac*scalefac/h_small - origin
        m1 = np.array(f[u'Header'].attrs[u'MassTable'][1]) * 10**10/h_small * np.ones(np.shape(xyz1)[0])
        dist1 = np.sqrt( xyz1[:,0]*xyz1[:,0] + xyz1[:,1]*xyz1[:,1] + xyz1[:,2]*xyz1[:,2] )

        m1 = np.concatenate(( m1, 1.e-3*np.ones_like(radArr), 1.e-3*np.ones_like(radArr) ))
        dist1 = np.concatenate(( dist1, radArr+1.e-6, radArr-1.e6 ))

    else:
        m1 = np.zeros_like(radArr)
        dist1 = np.zeros_like(radArr)


    if 4 in ptypes:

        xyz4 = np.array(f[u'PartType4'][u'Coordinates']) * dist_unit_fac*scalefac/h_small - origin
        m4 = np.array(f[u'PartType4'][u'Masses']) * 10**10/h_small
        dist4 = np.sqrt( xyz4[:,0]*xyz4[:,0] + xyz4[:,1]*xyz4[:,1] + xyz4[:,2]*xyz4[:,2] )

        m4 = np.concatenate(( m4, 1.e-3*np.ones_like(radArr), 1.e-3*np.ones_like(radArr) ))
        dist4 = np.concatenate(( dist4, radArr+1.e-6, radArr-1.e6 ))

    else:
        m4 = np.zeros_like(radArr)
        dist4 = np.zeros_like(radArr)

    mass_allpt = [ m0, m1, 'bigwhale', 'bigwhale', m4, 'bigwhale' ]
    dist_allpt = [ dist0, dist1, 'bigwhale', 'bigwhale', dist4, 'bigwhale' ]

    for ptype in ptypes:

        #ind = np.argsort(dist_allpt[ptype])
        #distOrd = dist_allpt[ptype][ind]
        #massOrd = mass_allpt[ptype][ind]
        ordered = order_masses( dist_allpt[ptype], mass_allpt[ptype] )
        distOrd = ordered[0]
        #cummass = ordered[1]

        indrads = np.array([ np.argmax(distOrd>rad_j) for rad_j in radArr ]) # index of 1st ptcl that exceeds measured radius
        distOrd = distOrd[indrads] # distOrd
        cummass = ordered[1][indrads] # cummass # sum masses of all existing particles first, then take the index
        #print('Should basically all be 1.e-6', distOrd-radArr)

        # calculate
        vol = 4.*np.pi/3. * distOrd*distOrd*distOrd
        density = cummass / vol
        velocity = np.sqrt( G * cummass / distOrd )

        ans[ptype,0] = distOrd
        ans[ptype,1] = cummass #cummass
        ans[ptype,2] = density #density
        ans[ptype,3] = velocity #velocity

        #print('All 4 must be the same:', len(radArr),len(cummass),len(density),len(velocity))

    ans[6,0] = radArr
    ans[6,1] = sum([ ans[i,1] for i in range(6) ]) #sum cumulative masses
    ans[6,2] = sum([ ans[i,2] for i in range(6) ]) #sum densities
    ans[6,3] = np.sqrt( sum([ ans[i,3]*ans[i,3] for i in range(6) ]) ) #sum velocities in quadrature

    return ans
