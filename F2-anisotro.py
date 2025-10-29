import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import zHead as hd
from astropy.cosmology import FlatLambdaCDM

arg_options = hd.handle_arguments()
ActiveRunIDs = arg_options[0]
print(ActiveRunIDs)
date = arg_options[1]
NRUNS = len(ActiveRunIDs)

print(hd.runlabels)

###

MODE = hd.DEFAULT_MODE
zetas = [1.0, 3.0, 10.0, 30.0, 0.3, 0.1, 1.0, 10.0, 1.0, 10.0]

###

gammas = []
gamma_errs = []
ptot_lins_all = []
ptot_vecs_all = []
ptot_lins_errs = []
ptot_vecs_errs = []
for j,(RunName, f) in enumerate(zip(hd.runlabels,hd.finalsnap_list)):
    if j in ActiveRunIDs:

        print(j,f)

        px = 1.e10 * np.array(f[u'PartType4'][u'CumInjSNMomentum_x'])
        py = 1.e10 * np.array(f[u'PartType4'][u'CumInjSNMomentum_y'])
        pz = 1.e10 * np.array(f[u'PartType4'][u'CumInjSNMomentum_z'])
        age = 2000. - 1000. * np.array(f[u'PartType4'][u'GFM_StellarFormationTime'])
        ptot_lin = px+py+pz
        ptot_vec = np.sqrt(px*px + py*py + pz*pz)

        ratio_zx = pz/px
        kk_notnan_zx = np.logical_not(np.isnan(ratio_zx))
        gamma = np.average( ratio_zx[kk_notnan_zx] )
        std_gamma = np.std( ratio_zx[kk_notnan_zx] )

        ratio_yx = py/px
        kk_notnan_yx = np.logical_not(np.isnan(ratio_yx))
        should_be_1 = np.average(ratio_yx[kk_notnan_yx])
        std_yx = np.std(ratio_yx[kk_notnan_yx])

        print( px )
        print( zetas[j], np.average(ptot_lin)/1.e7, np.average(ptot_vec)/1.e7, np.std(ptot_lin)/1.e7, np.std(ptot_vec)/1.e7 )
        print( zetas[j], gamma, std_gamma )

        gammas.append(gamma)
        gamma_errs.append(std_gamma)
        ptot_lins_all.append( np.average(ptot_lin) )
        ptot_vecs_all.append( np.average(ptot_vec) )
        ptot_lins_errs.append( np.std(ptot_lin) )
        ptot_vecs_errs.append( np.std(ptot_vec) )

    else:
        gammas.append(np.nan)
        gamma_errs.append(np.nan)

np.savetxt( 'data/gamma_dat', np.transpose(np.array([ zetas, gammas, gamma_errs, ptot_lins_all, ptot_vecs_all, ptot_lins_errs, ptot_vecs_errs ])) )

###

gammadat = np.loadtxt('data/gamma_dat')

### ANISOTROPY FIGURE 2

fig,ax = plt.subplots(1,1,figsize=(8.5,6.5))

zetaline = 10.**np.linspace(-2.5,1.75,50)

ax.plot( zetaline, zetaline, linestyle='dashed', linewidth=3.7, alpha=0.6, color='#458854', label='1:1' )

ax.fill_between( zetaline, 0.83*np.ones(50), y2=1.2*np.ones(50), alpha=0.37, color='#999999' )

ax.errorbar( zetas[:6], gammas[:6], yerr=gamma_errs[:6], color='black', fmt='o', markersize=12, elinewidth=2.5, ecolor='black' )

ax.text( 13.0, 1.0, 'Isotropic Injection', fontsize=16, color='#999999', ha='center', va='center' )

ax.text( 0.032, 1.5, r'$\uparrow$Out-of-Disk Favored', fontsize=16.5, color='#993355', ha='left', alpha=0.84 )
ax.text( 0.032, 0.67, r'$\downarrow$In-Disk Favored', fontsize=16.5, color='#553399', ha='left', va='top', alpha=0.84 )

ax.tick_params(which='both',direction='in',labelsize=16)
ax.legend(loc='lower right', fontsize=16)
ax.set_xlabel(r'Nominal Anisotropy $\zeta$', fontsize=17)
ax.set_ylabel(r'Effective Anisotropy $\gamma = \langle \Sigma p_z / \Sigma p_x \rangle_{*}$', fontsize=17)

ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')

ax.set_xlim(( 0.03, 50. ))
ax.set_ylim(( 0.05, 30. ))
ax.set_xscale('log')
ax.set_yscale('log')

fig.savefig('OfficialPaperPlots/F2-AnisotropyOfRuns.png')
plt.show()
