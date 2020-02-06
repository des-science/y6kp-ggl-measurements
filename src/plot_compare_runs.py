import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
from info import config, path_config
import twopoint


'''
Script to plot the comparison between two runs, i.e. compare two gammat runs.
'''

plot_gammat = 1
plot_ratios = 0
color_a = 'teal'
color_b = 'powderblue'

print '------------------------------------\nStarting plotting script...\n------------------------------------'
config['lens_v'] = 'redmagic'
config['thlims'] = np.array([2.5,250.])

config1 = dict(config)
config1['lens_w'] = True
path1 = path_config(config1)

config2 = dict(config)
config2['lens_w'] = False

path1 = path_config(config1)
path2 = path_config(config2)
print path1
print path2

def path_measurement(path):
    return '../runs/'+path+'/measurement/'

def path_plots(path):
    return '../plots/'+path+'/'


filename_a = 'gt_twopointfile_BLINDED'
filename_b = 'gt_boosted_twopointfile_BLINDED'

T1_a = twopoint.TwoPointFile.from_fits(path_measurement(path1)+ filename_a + '.fits') # no boosts, with weights
T2_a = twopoint.TwoPointFile.from_fits(path_measurement(path2)+ filename_a + '.fits') # no boosts, no weights
name_plot_a = path_plots(path1) + filename_a + '_lens_w_effect'

T1_b = twopoint.TwoPointFile.from_fits(path_measurement(path1)+ filename_b + '.fits') # boosts, with weights
T2_b = twopoint.TwoPointFile.from_fits(path_measurement(path2)+ filename_b + '.fits') # boosts, no weights
name_plot_b = path_plots(path1) + filename_b + '_lens_w_effect'


if plot_gammat:
    T1_a.plots(name_plot_a, plot_kernel = False, plot_cov = False, save_pickle = True, latex = False, label_legend = 'With LSS weights (no boosts)', blind_yaxis=True)
    T2_a.plots(name_plot_a, plot_kernel = False, plot_cov = False, load_pickle = True, latex = False, label_legend = 'No weights (no boosts)', blind_yaxis=True)

    T1_b.plots(name_plot_b, plot_kernel = False, plot_cov = False, save_pickle = True, latex = False, label_legend = 'With LSS weights (and boosts)', blind_yaxis=True)
    T2_b.plots(name_plot_b, plot_kernel = False, plot_cov = False, load_pickle = True, latex = False, label_legend = 'No weights (and boosts)', blind_yaxis=True)

if plot_ratios:

    ratio_errors = 0
    ratio_signals = 1
    
    s1_a = T1_a.spectra[0]
    s2_a = T2_a.spectra[0]

    s1_b = T1_b.spectra[0]
    s2_b = T2_b.spectra[0]

    pairs = s1_a.bin_pairs
    npairs = len(pairs)
    bins1 = np.transpose(pairs)[0]
    bins2 = np.transpose(pairs)[1]
    nbins1 = np.max(bins1)
    nbins2 = np.max(bins2)

    shade_until = None
    if ratio_errors:
        label = 'Err gt (With weights/No weights)'
    if ratio_signals:
        label = 'gt (With weights/No weights)'

    fig, ax = plt.subplots(nbins2, nbins1, figsize=(2.2*nbins1, 2.2*nbins2), sharey=True, sharex=True)
    for k,pair in enumerate(pairs):
            i,j = pair

            # No boosts
            theta1_a, xi1_a = s1_a.get_pair(i,j)
            error1_a = s1_a.get_error(i,j)
            theta2_a, xi2_a = s2_a.get_pair(i,j)
            error2_a = s2_a.get_error(i,j)
            assert theta1_a.all() ==theta2_a.all()

            if ratio_errors:
                ax[j-1][i-1].plot(theta1_a, error1_a/error2_a, color = color_a, label = 'No boosts')
            if ratio_signals:
                ax[j-1][i-1].plot(theta1_a, xi1_a/xi2_a, color = color_a, label = 'No boosts')

            # With boosts
            theta1_b, xi1_b = s1_b.get_pair(i,j)
            error1_b = s1_b.get_error(i,j)
            theta2_b, xi2_b = s2_b.get_pair(i,j)
            error2_b = s2_b.get_error(i,j)
            assert theta1_b.all() ==theta2_b.all()
            if ratio_errors:
                ax[j-1][i-1].plot(theta1_b, error1_b/error2_b, color = color_b, label = 'With boosts')
            if ratio_signals:
                ax[j-1][i-1].plot(theta1_b, xi1_b/xi2_b, color = color_b, label = 'With boosts')

            ax[j-1][i-1].axhline(y=1, color = 'k', ls=':')
            ax[j-1][i-1].text(0.85, 0.85, "{},{}".format(i, j), horizontalalignment='center',
                              verticalalignment='center', transform=ax[j-1][i-1].transAxes, fontsize=12)
            if shade_until is not None:
                if not load_pickle:
                    ax[j-1][i-1].axvspan(min(theta)*0.8, shade_until[i-1], color='gray', alpha=0.2)
                    ax[j-1][i-1].set_xlim(left=min(theta)*0.8, right=max(theta)*1.2)

            ax[j-1][i-1].set_xscale('log', nonposx='clip')
            ax[j-1][i-1].set_ylim(ymin=0.9, ymax=1.1)
            ax[j-1][i-1].xaxis.set_major_formatter(
                ticker.FormatStrFormatter('$%d$'))
            
            if (not all(bins1 == bins2)) & (j == nbins2):
                ax[j-1][i-1].set_xlabel(r"$\theta$ [arcmin]")
            if all(bins1 == bins2):
                ax[j-1][i-1].set_xlabel(r"$\theta$ [arcmin]")
            if i == 1:
                ax[j-1][i-1].set_ylabel(label, fontsize = 7)

    plt.legend(prop={'size': 8})
    if ratio_errors:
        name = name_plot_b + '_ratio_errors'
    if ratio_signals:
        name = name_plot_b + '_ratio_signals'
    plt.savefig(name + '.pdf')
    plt.savefig(name + '.png')

