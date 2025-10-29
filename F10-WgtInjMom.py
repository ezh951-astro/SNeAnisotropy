import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import zHead as hd

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

ptot_lins_all_pmf = []
ptot_vecs_all_pmf = []
ptot_lins_errs_pmf = []
ptot_vecs_errs_pmf = []

for j,(RunName, f) in enumerate(zip(hd.runlabels,hd.finalsnap_list)):
    if j in ActiveRunIDs:

        print(j,f)

        px = 1.e10 * np.array(f[u'PartType4'][u'CumInjSNMomentum_x'])
        py = 1.e10 * np.array(f[u'PartType4'][u'CumInjSNMomentum_y'])
        pz = 1.e10 * np.array(f[u'PartType4'][u'CumInjSNMomentum_z'])

        #age = 2000. - 1000. * np.array(f[u'PartType4'][u'GFM_StellarFormationTime'])

        ptot_lin = px+py+pz
        ptot_vec = np.sqrt(px*px + py*py + pz*pz)

        mformed = 1.e10 * np.array(f[u'PartType4'][u'GFM_InitialMass'])

        ptot_lin_permf = ptot_lin / mformed
        ptot_vec_permf = ptot_vec / mformed # units of km/s, probably on the order of 10^3-4 or so

        ratio_zx = pz/px
        kk_notnan_zx = np.logical_not(np.isnan(ratio_zx))
        gamma = np.average( ratio_zx[kk_notnan_zx] )
        std_gamma = np.std( ratio_zx[kk_notnan_zx] )

        #ratio_yx = py/px
        #kk_notnan_yx = np.logical_not(np.isnan(ratio_yx))
        #should_be_1 = np.average(ratio_yx[kk_notnan_yx])
        #std_yx = np.std(ratio_yx[kk_notnan_yx])

        #print( px )
        #print( zetas[j], np.average(ptot_lin)/1.e7, np.average(ptot_vec)/1.e7, np.std(ptot_lin)/1.e7, np.std(ptot_vec)/1.e7 )
        #print( zetas[j], gamma, std_gamma )

        gammas.append(gamma)
        gamma_errs.append(std_gamma)

        ptot_lins_all_pmf.append( np.average(ptot_lin_permf) )
        ptot_vecs_all_pmf.append( np.average(ptot_vec_permf) )
        ptot_lins_errs_pmf.append( np.std(ptot_lin_permf) )
        ptot_vecs_errs_pmf.append( np.std(ptot_vec_permf) )

gammadat =  np.transpose(np.array([ zetas, gammas, gamma_errs, ptot_lins_all_pmf, ptot_vecs_all_pmf, ptot_lins_errs_pmf, ptot_vecs_errs_pmf ]))
print(np.shape(gammadat))

### MOMENTUM INJECTION FIGURE FIG. 10

figM,axM = plt.subplots(2,1,figsize=(8.5,7.5))
plt.subplots_adjust(hspace=0.25)

axM[0].errorbar( gammas[:6], gammadat[:,3][:6], xerr=gammadat[:,2][:6], yerr=gammadat[:,5][:6], color='#AD1519', fmt='o', label='Solid-Angle Weight (Default)' )
axM[0].errorbar( gammas[6:8], gammadat[:,3][6:8], xerr=gammadat[:,2][6:8], yerr=gammadat[:,5][6:8], color='#004494', fmt='o', label='Volume Weight' )
axM[0].errorbar( gammas[8:], gammadat[:,3][8:], xerr=gammadat[:,2][8:], yerr=gammadat[:,5][8:], color='#006633', fmt='o', label='Mass Weight' )

axM[1].errorbar( gammas[:6], gammadat[:,4][:6], xerr=gammadat[:,2][:6], yerr=gammadat[:,6][:6], color='#AD1519', fmt='s', label='Solid-Angle Weight' )
axM[1].errorbar( gammas[6:8], gammadat[:,4][6:8], xerr=gammadat[:,2][6:8], yerr=gammadat[:,6][6:8], color='#004494', fmt='s', label='Volume Weight' )
axM[1].errorbar( gammas[8:], gammadat[:,4][8:], xerr=gammadat[:,2][8:], yerr=gammadat[:,6][8:], color='#006633', fmt='s', label='Mass Weight' )

axM[1].set_xlabel(r'Effective Anisotropy $\gamma$', fontsize=17)

axM[0].set_title('Linear', fontsize=16)
axM[1].set_title('Quadrature', fontsize=16)

figM.supylabel('Injected Momentum per $M_{*}$ Formed [km/s]', fontsize=17.5)

axM[0].set_ylim(( 6.e2, 1.5e4 ))
axM[1].set_ylim(( 6.e2, 1.5e4 ))

for k in range(2):
    axM[k].tick_params(labelsize=15, which='both', direction='in')
    axM[k].set_xlim(( 0.05, 20. ))
    axM[k].set_xscale('log')
    axM[k].set_yscale('log')

axM[1].legend(loc='lower left',fontsize=14.5)

#axM[0].set_xticklabels([])

axM[0].xaxis.set_ticks_position('both')
axM[1].xaxis.set_ticks_position('both')

axM[0].yaxis.set_ticks_position('both')
axM[1].yaxis.set_ticks_position('both')

figM.savefig('OfficialPaperPlots/F10-MomentumOfRuns.png')
plt.show()
