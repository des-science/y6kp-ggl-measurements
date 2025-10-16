import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from numpy.random import randint
from matplotlib import ticker
import twopoint

def radian_to_arcmin(theta_radian):
    return theta_radian*180/np.pi*60

output = {}
output['3x2_bf'] = {}
output['3x2_bf']['path'] = 'fiducial32pt_lin_lcdm_sim_bestfit_values/'
output['3x2_bf']['label'] = r'3x2pt best-fit' 

#data_file = '../../../y6-3x2pt/inference/data_vectors/simulated/sim_best_32pt_lcdm_3x2pt_2025-07-28-20h_UNBLINDED_wmask_nob2.fits'
data_file = '../../../y6-3x2pt/inference/data_vectors/real/3x2pt_2025-07-28-20h_UNBLINDED_wmask_nob2.fits'
#data_file = 'sim_y6.fits'

nbins_l = 6
nbins_s = 4

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fontsize = 16

######################################################################################
#definition of the plotting functions
##################################################################

def plot_cl_shear(dic, spectrum_names):
    ell = np.loadtxt(dic['path']+spectrum_name[0]+'ell.txt')
    for spectrum_name in spectrum_names:
        nbins1 = nbins2 = 4
        fig, ax = plt.subplots(nbins2, nbins1, figsize=(1.6*nbins1, 1.6*nbins2), sharey=False, sharex=True)
        for i in range(nbins1):
            for j in range(nbins2):
                if j > i:
                    fig.delaxes(ax[i, j])

                else:
                    print(i, j)
                    cl = np.loadtxt(dic['path']+spectrum_name+'bin_%d_%d.txt'%(i+1,j+1))
                    ax[i][j].plot(ell, cl, label = spectrum_name[0:-1])
                    ax[i][j].set_xlim(10,10000)
                    ax[i][j].set_xscale('log')

                    if i == (nbins1-1):
                        ax[i][j].set_xlabel('$\ell$')

    plt.title(dic['label'])
    plt.legend(loc=0)
    plt.show()
    plt.savefig('plots/prova.png')#output_dir+spectrum_name[0:-1]+'.png',format='png')
    plt.close()
    
def plot_cl_gglensing(dic, spectrum_names, name):

    
    #colors = ['k', 'k', 'teal', 'orange', 'gold', 'powderblue']
    colors = ['k',  'teal', 'orange', 'powderblue']
    ls = ['-',  '-', '-', '--']
    lw = [2., 2., 2., 2., 1.5, 1.5]
    fig, ax = plt.subplots(nbins_s, nbins_l, figsize=(3.1*nbins_l, 3.1*nbins_s), sharey=False, sharex=True)

    ell = np.loadtxt(dic['path']+spectrum_names[0]+'/ell.txt')
    for s, spectrum_name in enumerate(spectrum_names):
        for i in range(nbins_l):
            for j in range(nbins_s):
                cl = np.loadtxt(dic['path']+spectrum_name+'/bin_%d_%d.txt'%(i+1,j+1))
                ax[j][i].plot(ell, cl, color = colors[s], ls=ls[s], lw=lw[s], label = terms_ggl_cl_labels[s])
                ax[j][i].set_xlim(10,10000)
                ax[j][i].set_xscale('log')
                ax[j][i].axhline(y=0, ls = ':', color = 'k')
                #ax[j][i].set_yscale('log', nonposy='clip')
                ax[j][i].text(0.85, 0.85, "{},{}".format(i+1, j+1), horizontalalignment='center',
                                  verticalalignment='center', transform=ax[j][i].transAxes, fontsize=fontsize)
                
                ax[j][i].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                if i == 0:
                    ax[j][i].set_ylabel('$\mathcal{C}(l)$', fontsize=fontsize)
                if j == (nbins_s-1):
                    ax[j][i].set_xlabel('$l$', fontsize =fontsize)
                ax[j][i].tick_params(axis='both', which='major', labelsize=13)

    #ticker.ScalarFormatter(useOffset=True)
    legend = ax[0][2].legend(loc='upper center', bbox_to_anchor=(0.9, 1.4),
                             frameon= False, ncol=4, fontsize = fontsize)
    #legend = ax[0][2].legend(loc='upper center', bbox_to_anchor=(0.5, 1.4),
    #          frameon= False, ncol=4, title=dic['label'], fontsize = 13)
    plt.setp(legend.get_title(),fontsize='13')

    plt.savefig('plots/gglensing_terms_fourier_%s.png'%name, dpi = 500, bbox_inches = 'tight',  pad_inches = 0.1)
    plt.savefig('plots/gglensing_terms_fourier_%s.pdf'%name, bbox_inches = 'tight',  pad_inches = 0.1)
    plt.close()
    
def plot_data(ax):
    #T1 = twopoint.TwoPointFile.from_fits('/Users/juditprat/Downloads/2pt_NG_final_2ptunblind_10_26_20_wnz.fits')
    T1 = twopoint.TwoPointFile.from_fits(data_file)
    spectra = T1.spectra[2]
    
    for i in range(nbins_l):
        for j in range(nbins_s):
            theta1, xi1 = spectra.get_pair(i+1,j+1)
            error1 = spectra.get_error(i+1,j+1)
            print(i, j, len(theta1), len(xi1), len(error1))
            ax[j][i].errorbar(theta1, xi1, yerr=error1, fmt = 'o', color = 'k', alpha = 0.7, markersize = 3., capsize=1.5, capthick=0.8, elinewidth=0.8)
    return ax
    
def plot_xi_gglensing(dic, spectrum_names, name):
    colors = ['k', 'teal', 'orange', 'powderblue']
    ls = ['-', '-', '-', '--']
    lw = [2., 2., 2., 2., 1.5, 1.5]
    theta_cut = [26.83313651, 250, 13.61215672, 11.32891161, 10.01217238, 9.] # arcmin, corresponding to 6 Mpc/hw
    #theta_cut = [26.83313651, 17.63634989, 13.61215672, 11.32891161, 250, 250] # arcmin, corresponding to 6 Mpc/hw

    fig, ax = plt.subplots(nbins_s, nbins_l, figsize=(3.1*nbins_l, 3.1*nbins_s), sharey=False, sharex=True)
    theta = radian_to_arcmin(np.loadtxt(dic['path']+spectrum_names[0]+'/theta.txt'))
    mask = (theta<250)&(theta>0.25)

    ax = plot_data(ax)
    
    for s, spectrum_name in enumerate(spectrum_names):
        for i in range(nbins_l):
            for j in range(nbins_s):
                xi = np.loadtxt(dic['path']+spectrum_name+'/bin_%d_%d.txt'%(i+1,j+1))
                ax[j][i].plot(theta[mask], xi[mask], color = colors[s], ls=ls[s], lw=lw[s], label = terms_ggl_xi_labels[s])
                ax[j][i].set_xlim(2.5,250)
                ax[j][i].set_xscale('log')
                #ax[j][i].set_yscale('log', nonposy='clip')
                ax[j][i].text(0.85, 0.85, "{},{}".format(i+1, j+1), horizontalalignment='center',
                                  verticalalignment='center', transform=ax[j][i].transAxes, fontsize=fontsize)

                ax[j][i].axvspan(2.5, theta_cut[i], color='lightgray', alpha=0.2)
                ax[j][i].axhline(y=0, ls = ':', color = 'k')
                ax[j][i].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                if i == 0:
                    ax[j][i].set_ylabel(r'$\gamma_t(\theta)$ and terms', fontsize = fontsize)
                if j == (nbins_s-1):
                    ax[j][i].set_xlabel(r'$\theta$ [arcmin]', fontsize = fontsize)
                ax[j][i].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%d$'))
                ax[j][i].tick_params(axis='both', which='major', labelsize=13)

    legend = ax[0][2].legend(loc='upper center', bbox_to_anchor=(0.9, 1.4),
                             frameon= False, ncol=4, fontsize = fontsize)
    #legend = ax[0][2].legend(loc='upper center', bbox_to_anchor=(0.5, 1.6),
    #          frameon= False, ncol=3, title=dic['label'], fontsize = 13)
    plt.setp(legend.get_title(),fontsize='16')
    plt.savefig('plots/gglensing_terms_real_data_%s.png'%name, dpi = 500, bbox_inches = 'tight',  pad_inches = 0.1)
    plt.savefig('plots/gglensing_terms_real_data_%s.pdf'%name, bbox_inches = 'tight',  pad_inches = 0.1)
    plt.close()
    


def plot_xi_gglensing_compare(dic, spectrum_names, name1, name2):
    colors = ['k', 'teal', 'orange', 'powderblue']
    lw = [2., 2., 2., 2., 1.5, 1.5]
    fig, ax = plt.subplots(nbins_s, nbins_l, figsize=(2.4*nbins_l, 2.4*nbins_s), sharey=False, sharex=True)
    theta1 = radian_to_arcmin(np.loadtxt(dic[name1]['path']+spectrum_names[0]+'/theta.txt'))
    theta2 = radian_to_arcmin(np.loadtxt(dic[name2]['path']+spectrum_names[0]+'/theta.txt'))
    assert (theta1==theta2).all()
    mask = (theta1<250)&(theta1>2.5)
    for s, spectrum_name in enumerate(spectrum_names):
        for i in range(nbins_l):
            for j in range(nbins_s):
                print(i, j)
                xi1 = np.loadtxt(dic[name1]['path']+spectrum_name+'/bin_%d_%d.txt'%(i+1,j+1))
                xi2 = np.loadtxt(dic[name2]['path']+spectrum_name+'/bin_%d_%d.txt'%(i+1,j+1))
                ax[j][i].plot(theta1[mask], xi1[mask], color = colors[s], ls='-', lw=lw[s], label = spectrum_name + ' ' + name1)
                ax[j][i].plot(theta2[mask], xi2[mask], color = colors[s], ls='--', lw=lw[s], label = spectrum_name + ' ' + name2)
                ax[j][i].set_xlim(2.5,250)
                ax[j][i].set_xscale('log')
                ax[j][i].text(0.85, 0.85, "{},{}".format(i+1, j+1), horizontalalignment='center',
                                  verticalalignment='center', transform=ax[j][i].transAxes, fontsize=12)

                ax[j][i].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                if i == 0:
                    ax[j][i].set_ylabel(r'$\gamma_t(\theta)$ and terms')
                if j == (nbins_s-1):
                    ax[j][i].set_xlabel(r'$\theta$ [arcmin]')

    #ticker.ScalarFormatter(useOffset=True)
    #fig.suptitle(dic['label'])
    legend = ax[0][2].legend(loc='upper center', bbox_to_anchor=(0.9, 1.6),
                             frameon= False, ncol=3, title=dic[name1]['label'] + ',' + dic[name2]['label'], fontsize = 16)
    plt.setp(legend.get_title(),fontsize='13')

    plt.savefig('plots/gglensing_terms_real_compare_%s_%s.png'%(name1, name2), dpi=400)
    plt.close()
    



#############################
# Lucas code below

def plot_cl_ratio(spectrum_name):
    bins = listdir(dv1_dir+spectrum_name)
    bins.remove('ell.txt')
    bins.remove('values.txt')
    dv1_ell = np.loadtxt(dv1_dir+spectrum_name+'ell.txt')
    dv2_ell = np.loadtxt(dv2_dir+spectrum_name+'ell.txt')
    print( 'is dv1_ell == dv2_ell?',all(dv1_ell==dv2_ell))
    pl.figure()
    for bin_name in bins:
        dv1 = np.loadtxt(dv1_dir+spectrum_name+bin_name)
        dv2 = np.loadtxt(dv2_dir+spectrum_name+bin_name)
        pl.plot(dv1_ell, dv1/dv2, ls=linestyles[randint(6)], color=colors[randint(7)],label=bin_name)
        
    pl.xlim(10,10000)
    pl.ylim(0.95,1.05)
    pl.xscale('log')
    pl.axhline(y=1,color='k')
    pl.xlabel('ell')
    pl.title('dv1/dv2 for '+spectrum_name)
    pl.legend(loc=0)
    #pl.show()
    pl.savefig(output_dir+spectrum_name[0:-1]+'.png',format='png')
    pl.close()

def plot_power(spectrum_name): #spectrum_name are intrinsic_power, matter_power_nl, matter_power_lin, matter_intrinsic_power....
    dv2_k = np.loadtxt(dv2_dir+spectrum_name+'k_h.txt')
    dv2_pk = np.loadtxt(dv2_dir+spectrum_name+'p_k.txt')
    dv1_k = np.loadtxt(dv1_dir+spectrum_name+'k_h.txt')
    dv1_pk = np.loadtxt(dv1_dir+spectrum_name+'p_k.txt')

    print('is dv1_k == dv2_k?',all(dv1_k==dv2_k))
    pl.figure()
    pl.plot(dv2_k, dv1_pk[0]/dv2_pk[0], ls=linestyles[randint(6)], color=colors[randint(7)],label='redshift = 0')
    pl.xscale('log')
    pl.xlim(0.0001,10)
    pl.ylim(0.95,1.05)
    pl.axhline(y=1,color='k')
    pl.xlabel('k')
    pl.title('dv1/dv2 for '+spectrum_name)
    pl.legend(loc=0)
    #pl.show()
    pl.savefig(output_dir+spectrum_name[0:-1]+'.png',format='png')
    pl.close()

def plot_xi(spectrum_name):
    bins = listdir(dv2_dir+spectrum_name)
    bins.remove('theta.txt')
    bins.remove('values.txt')
    dv2_theta = np.loadtxt(dv2_dir+spectrum_name+'theta.txt')*(180.0*60.0/np.pi) #from radians to arcmin
    dv1_theta = np.loadtxt(dv1_dir+spectrum_name+'theta.txt')*(180.0*60.0/np.pi)

    print('is dv1_theta == dv2_theta?',all(dv1_theta==dv2_theta))
    
    pl.figure()
    for bin_name in bins:
        dv1 = np.loadtxt(dv1_dir+spectrum_name+bin_name)
        dv2 = np.loadtxt(dv2_dir+spectrum_name+bin_name)
        pl.plot(dv2_theta, dv1/dv2, ls=linestyles[randint(6)], color=colors[randint(7)],label=bin_name)

    
    pl.xlim(2.5,250)
    pl.ylim(0.95,1.05)
    pl.xscale('log')
    pl.axhline(y=1,color='k')
    pl.xlabel('theta [arcmin]')
    pl.title('dv1/dv2 for '+spectrum_name)
    pl.legend(loc=0)
    #pl.show()
    pl.savefig(output_dir+spectrum_name[0:-1]+'.png',format='png')
    pl.close()
    
###################################################################


terms_clustering = ['galaxy_cl', 'galaxy_cl_gg', 'magnification_cl', 'galaxy_magnification_cl'] 
terms_ggl_cl = ['galaxy_shear_cl', 'magnification_shear_cl', 'galaxy_intrinsic_cl', 'magnification_intrinsic_cl' ]
terms_ggl_cl_labels = [r'$C_{gm}$', r'$2(\alpha-1)\,  C_{mm}$ (Lens mag.)', r'$C_{gI}$ (IA)', r'$C_{mI}$ (Lens mag. $\times$ IA)']
terms_ggl_xi = ['galaxy_shear_xi', 'magnification_shear_xi', 'galaxy_intrinsic_xi', 'magnification_intrinsic_xi' ]
terms_ggl_xi_labels = [r'Total $\gamma_t$', 'Lens magnification', 'Intrinsic Alignments', r'Lens magnification $\times$ Intrinsic Alignments']
terms_shear = ['shear_cl']


plot_cl_gglensing(output['3x2_bf'], terms_ggl_cl, '3x2_bf')
plot_xi_gglensing(output['3x2_bf'], terms_ggl_xi, '3x2_bf')

#plot_xi_gglensing_compare(output, terms_ggl_xi, 'NLA', 'TATT5x')
#plot_cl_gglensing(output['TATT_onlyA2neg'], terms_ggl_cl, 'TATT_onlyA2neg')
#plot_xi_gglensing(output['TATT_onlyA2neg'], terms_ggl_xi, 'TATT_onlyA2neg')



'''
#plot_cl('shear_cl_ii/')
#plot_cl('shear_cl_gi/')
plot_cl('galaxy_shear_cl/')
plot_cl('galaxy_intrinsic_cl/')
plot_cl('galaxy_cl/')

plot_power('matter_power_lin/')
plot_power('matter_power_nl/')
plot_power('matter_galaxy_power/')
plot_power('intrinsic_power/')
plot_power('matter_intrinsic_power/')
plot_power('galaxy_power/')
plot_power('galaxy_intrinsic_power/')

#plot_xi('galaxy_intrinsic_xi/')
plot_xi('galaxy_shear_xi/')
plot_xi('galaxy_xi/')
plot_xi('shear_xi_plus/')
plot_xi('shear_xi_minus/')
#plot_xi('shear_xi_gi/')
#plot_xi('shear_xi_ii/')
'''
