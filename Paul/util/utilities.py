## MAKE SURE YOU SET THE ENVIRONMENT VARIABLE 'PYDIR' = directory with all the PFH routines
##
## this loads various parts of my python libraries under one 'parent' module heading
##  -- this is purely for convenience, any of these can be called on their own --
##

import numpy as np
import math
import colors_sps.colors_table as ctab
import colors_sps.lum_mag_conversions as lmag
import attenuation.attenuate_wrapper as atten
import attenuation.cross_section as atten_cx
import agn_spectrum.agn_spectrum_wrapper as agnspec

def return_python_routines_homedir():
    import os
    return os.environ['PYDIR']
    ## make sure this is set to the pfh_python routines directory, 
    ##   as an environment variable, since everything will link through it! otherwise nothing works
    ## alternatively, you can manually set this here, but you must remember to change the code 
    ##   on every machine you install onto !
    ##return '/Users/phopkins/Documents/work/plots/python'

def return_python_routines_cdir():
    return return_python_routines_homedir()+'/c_libraries'
    
## routines from colors_sps module
def colors_table( age_in_Gyr, metallicity_in_solar_units, 
    BAND_ID=0, SALPETER_IMF=0, CHABRIER_IMF=1, QUIET=0, CRUDE=0, 
    RETURN_NU_EFF=0, RETURN_LAMBDA_EFF=0, UNITS_SOLAR_IN_BAND=0 ):
    return ctab.colors_table( age_in_Gyr, metallicity_in_solar_units, 
        BAND_ID=BAND_ID, SALPETER_IMF=SALPETER_IMF, CHABRIER_IMF=CHABRIER_IMF, QUIET=QUIET, CRUDE=CRUDE, 
        RETURN_NU_EFF=RETURN_NU_EFF, RETURN_LAMBDA_EFF=RETURN_LAMBDA_EFF, UNITS_SOLAR_IN_BAND=UNITS_SOLAR_IN_BAND )

def get_solar_mags():
    return lmag.get_solar_mags();
    
def luminosity_to_magnitude( L, \
	UNITS_SOLAR_BOL=0, UNITS_SOLAR_BAND=0, \
	UNITS_CGS=0, \
	NU_L_NU=0, L_NU=0, \
	BAND_U=0,BAND_B=0,BAND_V=0,BAND_R=0,BAND_I=0, \
	BAND_J=0,BAND_H=0,BAND_K=0,BAND_SDSS_u=0, \
	BAND_SDSS_g=0,BAND_SDSS_r=0,BAND_SDSS_i=0, \
	BAND_SDSS_z=0,BOLOMETRIC=0, BAND_NUMBER=0, \
	VEGA=0, AB=0 ):
    return lmag.luminosity_to_magnitude( L, 
	UNITS_SOLAR_BOL=UNITS_SOLAR_BOL, UNITS_SOLAR_BAND=UNITS_SOLAR_BAND, 
	UNITS_CGS=UNITS_CGS, 
	NU_L_NU=NU_L_NU, L_NU=L_NU, 
	BAND_U=BAND_U,BAND_B=BAND_B,BAND_V=BAND_V,BAND_R=BAND_R,BAND_I=BAND_I,
	BAND_J=BAND_J,BAND_H=BAND_H,BAND_K=BAND_K,BAND_SDSS_u=BAND_SDSS_u, 
	BAND_SDSS_g=BAND_SDSS_g,BAND_SDSS_r=BAND_SDSS_r,BAND_SDSS_i=BAND_SDSS_i, 
	BAND_SDSS_z=BAND_SDSS_z,BOLOMETRIC=BOLOMETRIC, BAND_NUMBER=BAND_NUMBER, 
	VEGA=VEGA, AB=AB, MAGNITUDE_TO_LUMINOSITY=0)
	
def magnitude_to_luminosity( L, \
	UNITS_SOLAR_BOL=0, UNITS_SOLAR_BAND=0, \
	UNITS_CGS=0, \
	NU_L_NU=0, L_NU=0, \
	BAND_U=0,BAND_B=0,BAND_V=0,BAND_R=0,BAND_I=0, \
	BAND_J=0,BAND_H=0,BAND_K=0,BAND_SDSS_u=0, \
	BAND_SDSS_g=0,BAND_SDSS_r=0,BAND_SDSS_i=0, \
	BAND_SDSS_z=0,BOLOMETRIC=0, BAND_NUMBER=0, \
	VEGA=0, AB=0 ):
    return lmag.luminosity_to_magnitude( L, 
	UNITS_SOLAR_BOL=UNITS_SOLAR_BOL, UNITS_SOLAR_BAND=UNITS_SOLAR_BAND, 
	UNITS_CGS=UNITS_CGS, 
	NU_L_NU=NU_L_NU, L_NU=L_NU, 
	BAND_U=BAND_U,BAND_B=BAND_B,BAND_V=BAND_V,BAND_R=BAND_R,BAND_I=BAND_I,
	BAND_J=BAND_J,BAND_H=BAND_H,BAND_K=BAND_K,BAND_SDSS_u=BAND_SDSS_u, 
	BAND_SDSS_g=BAND_SDSS_g,BAND_SDSS_r=BAND_SDSS_r,BAND_SDSS_i=BAND_SDSS_i, 
	BAND_SDSS_z=BAND_SDSS_z,BOLOMETRIC=BOLOMETRIC, BAND_NUMBER=BAND_NUMBER, 
	VEGA=VEGA, AB=AB, MAGNITUDE_TO_LUMINOSITY=1)
	
	
## routines from attenuation module
def attenuate( nu_in_Hz, log_NH, metallicity_in_solar, \
        SMC=0, LMC=0, MW=0, BB=0, IR=0, SX=0, HX=0):
    return atten.attenuate( nu_in_Hz, log_NH, metallicity_in_solar, \
        SMC=SMC, LMC=LMC, MW=MW, BB=BB, IR=IR, SX=SX, HX=HX)
        
def cross_section( f0, METALLICITY_OVER_SOLAR=1.0 ):
    return atten_cx.cross_section( f0, METALLICITY_OVER_SOLAR=METALLICITY_OVER_SOLAR);
       
def opacity_per_solar_metallicity( f0 ):
    return atten_cx.opacity_per_solar_metallicity( f0 );
       
## routines from agn_spectrum module
def agn_spectrum( nu_in_Hz, log_l_bol, \
	BB=0, IR=0, SX=0, HX=0, \
	HRH=0, MARCONI=0, RICHARDS=0, \
	SDSS=0 ):
    return agnspec.agn_spectrum( nu_in_Hz, log_l_bol, \
	BB=BB, IR=IR, SX=SX, HX=HX, \
	HRH=HRH, MARCONI=MARCONI, RICHARDS=RICHARDS, \
	SDSS=SDSS )


## procedure which will return, for a given input vector A_in, 
##    the perpendicular unit vectors B_out and C_out 
##    which form perpendicular axes to A
def return_perp_vectors(a, LOUD=0):
    eps = 1.0e-10
    a = np.array(a,dtype='f');
    a /= np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    for i in range(len(a)):
        if (a[i]==0.): a[i]=eps;
        if (a[i]>=1.): a[i]=1.-eps;
        if (a[i]<=-1.): a[i]=-1.+eps;
    a /= np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    ax=a[0]; ay=a[1]; az=a[2];

    ## use a fixed rotation of the a-vector by 90 degrees:
    ## (this anchors the solution so it changes *continously* as a changes)
    t0=np.double(math.pi/2.e0);
    bx=0.*ax; by=np.cos(t0)*ay-np.sin(t0)*az; bz=np.sin(t0)*ay+np.cos(t0)*az;
    ## c-sign is degenerate even for 'well-chosen' a and b: gaurantee right-hand 
    ##   rule is obeyed by defining c as the cross product: a x b = c
    cx=(ay*bz-az*by); cy=-(ax*bz-az*bx); cz=(ax*by-ay*bx); 
    B_out=np.zeros(3); C_out=np.zeros(3);
    B_out[:]=[bx,by,bz]; C_out[:]=[cx,cy,cz];
    
    if (LOUD==1):
        print( a  )
        print( B_out )  
        print( C_out )
        print( 'a_tot=',ax*ax+ay*ay+az*az) 
        print( 'b_tot=',bx*bx+by*by+bz*bz)
        print( 'c_tot=',cx*cx+cy*cy+cz*cz)
        print( 'adotb=',ax*bx+ay*by+az*bz)
        print( 'adotc=',ax*cx+ay*cy+az*cz)
        print( 'bdotc=',bx*cx+by*cy+bz*cz)
    
    return B_out, C_out
    
## vector nan checker (because I like dividing by zero, apparently)
def isnan(x):
    return np.isnan(x);

## more robust length checker (will catch scalars)
def checklen(x):
    return len(np.array(x,ndmin=1));

## round which guarantees an int result
def int_round(x):
    return np.int(np.round(x));

## generic tool for checking if values are ok
def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        return (np.isnan(input)==False) & (abs(input)<=xmax) & (input > 0.);
    else:
        return (np.isnan(input)==False) & (abs(input)<=xmax);

## return age of universe (for a flat universe) to a given redshift
def age_of_universe(z,h=0.71,Omega_M=0.27):
    ## exact solution for a flat universe
    a=1./(1.+z); x=Omega_M/(1.-Omega_M) / (a*a*a);
    t=(2./(3.*np.sqrt(1.-Omega_M))) * np.log( np.sqrt(x) / (-1. + np.sqrt(1.+x)) );
    t *= 13.777 * (0.71/h); ## in Gyr
    return t;

## return lookback time (for a flat universe) to a given redshift
def lookback_time(z,h=0.71,Omega_M=0.27):
    return age_of_universe(0.,h=h,Omega_M=Omega_M)-age_of_universe(z,h=h,Omega_M=Omega_M);

## simple moving-window signal smoothing (from scipy cookbook)
def smooth(x,window_len=11,window='hanning'):
        if x.ndim != 1:
                print("smooth only accepts 1 dimension arrays." )
                sys.exit()
        if x.size < window_len:
                print("Input vector needs to be bigger than window size." )
                sys.exit()
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                print("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'" )
                sys.exit()
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=np.ones(window_len,'d')
        else:  
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]

## linear interpolation, but with linear extrapolation for points beyond the data boundaries
##  (x_old is monotonically increasing, here)
def interp_w_extrap(x_new, x_old, y_old):
    y = np.interp(x_new,x_old,y_old)
    N = x_old.size - 1
    out = (x_new < x_old[0])
    if (y[out].size > 0):
        y[out] = y_old[0] + (y_old[1]-y_old[0])/(x_old[1]-x_old[0]) * (x_new[out]-x_old[0])
    out = (x_new > x_old[N])
    if (y[out].size > 0):
        y[out] = y_old[N] + (y_old[N]-y_old[N-1])/(x_old[N]-x_old[N-1]) * (x_new[out]-x_old[N])
    return y
    
## simple derivative function for y = y(x) [ x monotonically increasing ]
def my_simple_derivative(y_in, x_in):
    s = np.argsort(x_in)
    x = x_in[s]
    y = y_in[s]
    xm=x[0:x.size-1]
    xp=x[1:x.size]
    x_mid = 0.5*(xm+xp)
    dy_dx = np.diff(y)/np.diff(x)
    dy_dx_xin = interp_w_extrap( x, x_mid, dy_dx )
    z = 0.*y_in
    z[s] = dy_dx_xin
    return z
    
## return color based on an input number
def set_color(i):
    colors_vec=['black','blue','red','green','deepskyblue',\
        'darkviolet','orange','magenta','gold','sienna','pink',\
        'forestgreen','darkgray']
    if(i<0): i=0;
    if(i>=np.array(colors_vec).size): i=-1;
    return colors_vec[i];
    
## return marker based on an input number
def set_symbol(i):
    symbol_vec=['o','s','p','^','*','x','D','+',\
        'h','1','v','2','<','3','>','4','H']
    if(i<0): i=0;
    if(i>=np.array(symbol_vec).size): i=-1;
    return symbol_vec[i];

# data object contains data.pos, data.vel, and data.mass
# two key keywords:  edge_on and face_on
def determine_rotation_angle(data, edge_on=False, face_on=False, am_weight='sfr', **kwargs):
    if (kwargs.get("phi") is not None) or (kwargs.get("theta") is not None):
        phi = 0
        theta = 0
        if (kwargs.get("phi") is not None):   phi   = kwargs.get("phi")
        if (kwargs.get("theta") is not None): theta = kwargs.get("theta")
        return phi, theta

    if edge_on or face_on:
        x = data.pos[:,0]
        y = data.pos[:,1]
        z = data.pos[:,2]
        vx = data.vel[:,0]
        vy = data.vel[:,1]
        vz = data.vel[:,2]
        m = data.mass
        try: am_weight = getattr( data, am_weight )
        except: am_weight= np.ones_like(x)

        r = np.sqrt( x**2 + y**2 + z**2 )
        index = (am_weight > 0) & (r < 5.0)
        lz = np.sum( m[index] * (x[index] * vy[index] - y[index] * vx[index] ) )
        lx = np.sum( m[index] * (y[index] * vz[index] - z[index] * vy[index] ) )
        ly = np.sum( m[index] * (z[index] * vx[index] - x[index] * vz[index] ) )

        if face_on:
            phi   = np.arctan2( ly, lx )    #  + 3.14159/2.0
            theta =  np.arctan2( np.sqrt(lx**2 + ly**2), lz )      # + 3.14159/2.0
        if edge_on:
            phi   = np.arctan2( ly, lx ) + 3.14159/2.0
            theta = 3.14159/2.0 + np.arctan2( np.sqrt(lx**2 + ly**2), lz )      # + 3.14159/2.0
    else:
        phi = 0
        theta = 0
    return phi, theta

def rotate_data( data, phi=0, theta=0, **kwargs):
        x = data.pos[:,0]
        y = data.pos[:,1]
        z = data.pos[:,2]
        vx = data.vel[:,0]
        vy = data.vel[:,1]
        vz = data.vel[:,2]

        x_ = -z  * np.sin(theta) + (x * np.cos(phi) + y *np.sin(phi)) * np.cos(theta)
        y_ = -x  * np.sin(phi)   + y  * np.cos(phi)
        z_ =  z  * np.cos(theta) + (x * np.cos(phi) + y *np.sin(phi)) * np.sin(theta)
        vx_ = -vz  * np.sin(theta) + (vx * np.cos(phi) + vy *np.sin(phi)) * np.cos(theta)
        vy_ = -vx  * np.sin(phi)   + vy  * np.cos(phi)
        vz_ =  vz  * np.cos(theta) + (vx * np.cos(phi) + vy *np.sin(phi)) * np.sin(theta)

        data.pos[:,0] = x_
        data.pos[:,1] = y_
        data.pos[:,2] = z_

        data.vel[:,0] = vx_
        data.vel[:,1] = vy_
        data.vel[:,2] = vz_

        return data


def determine_image_stretch( data, set_dynrng=None, set_maxden=None, **kwargs):
    data.set_dynrng = set_dynrng
    data.set_maxden = set_maxden

    return data

def determine_image_bounds( data, xrange=None, yrange=None, zrange=None, **kwargs):
    if ( xrange is not None):               # assume here that we use manually prescribed values
        if (yrange is None): yrange= xrange         # if yrange is set, use that value.  Otherwise, copy xrange
        if (zrange is None): zrange = xrange        # same for zr
    else:
        print("No bounds detected [xr,yr,zr] in 'utilities.determine_image_bounds'.  Using min/max values" )
        x_stretch = np.max( data.pos[:,0]) - np.min(data.pos[:,0])
        y_stretch = np.max( data.pos[:,1]) - np.min(data.pos[:,1])
        x_mid     = 0.5*( np.max( data.pos[:,0]) + np.min(data.pos[:,0]) )
        y_mid     = 0.5*( np.max( data.pos[:,1]) + np.min(data.pos[:,1]) )
        z_mid     = 0.5*( np.max( data.pos[:,2]) + np.min(data.pos[:,2]) )

        stretch = np.max( [x_stretch, y_stretch] )

        xrange = [ x_mid - stretch/2.0 , x_mid + stretch/2.0 ]
        yrange = [ y_mid - stretch/2.0 , y_mid + stretch/2.0 ]
        zrange = [ z_mid - stretch/2.0 , z_mid + stretch/2.0 ]


    data.xr = xrange
    data.yr = yrange
    data.zr = zrange
    return data

def clip_particles( data, xr, yr, zr, **kwargs ):
    print ("Clipping paticles outside of range xr=[{:8.4f},{:8.4f}] and yr=[{:8.4f},{:8.4f}] and zr=[{:8.4f},{:8.4f}]".format( xr[0], xr[1], yr[0], yr[1], zr[0], zr[1] ) ) 

    okay_index = ( data.pos[:,0] > xr[0] ) & \
                 ( data.pos[:,0] < xr[1] ) & \
                 ( data.pos[:,1] > yr[0] ) & \
                 ( data.pos[:,1] < yr[1] ) & \
                 ( data.pos[:,2] > zr[0] ) & \
                 ( data.pos[:,2] < zr[1] ) 


    for field in dir(data):
        if type(getattr(data,field)) is np.ndarray:
            tmp = getattr(data,field)
            setattr( data, field, tmp[okay_index] ) 

    return data


