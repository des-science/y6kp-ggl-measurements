import numpy as np
import pyfits as pf
from info import zbins
import matplotlib.pyplot as plt

def get_zbins_match_nz(zs, nz_zlims, nzs):
    """
    To match the redshift distribution of a catalog to another. 
    zs: redshifts of the catalog you want to get matched. Eg. here MICE true redshifts.
    nz_zims: z array with form [zhist_low0, zhist_low1, ..., zhist_low-1, zhist_high-1].
    nzs: The nzs you want to match the distributions to.
    """
    print 'Starting function.' 
    zbins = -1*np.ones_like(zs, dtype=int)
    nz_weights = np.zeros_like(zs, dtype=float)
    bin_nums=np.arange(len(nzs))
    print 'Starting 1st loop.'
    for i,(zlo,zhi) in enumerate(zip(nz_zlims[:-1], nz_zlims[1:])):
        print i
        use = (zs>zlo)*(zs<=zhi)
        if use.sum()==0:
            continue
        draw_lims=nzs[:i]/nzs[:i].sum()
        zbins[use]=np.random.choice(bin_nums, size=use.sum(), p=nzs[:,i]/np.sum(nzs[:,i]))

    print 'Binning is done. Starting 2nd loop.'
    #now we've done the binning, let's also assign weights such that the n(z)s match.
    for zbin_i in bin_nums:
        print zbin_i
        use = zbins==zbin_i
        nz,_ = np.histogram(zs[use], bins=nz_zlims)
        nz_weights_by_bin = nzs[zbin_i].astype(float)/nz
        bin_inds=np.digitize(zs[use], nz_zlims)-1
        nz_weights[use] = nz_weights_by_bin[bin_inds]

    print "min/max weight:", np.min(nz_weights), np.max(nz_weights)
    return zbins, nz_weights

# Load MICE
path_mice = '../../../y1_shear_tests/sims/mice/'
shear = pf.getdata(path_mice + 'mice2_shear.fits')
ids = shear['unique_gal_id']
ra = shear['ra_gal']
dec = shear['dec_gal']
z_true = shear['z_cgal']
g1 = shear['gamma1']
g2 = shear['gamma2']

# If downsample is true it select only 1% of catalog randomly
downsample = False
if downsample:
    mask_random = np.random.choice([False, True], len(ra), p=[0.99, 0.01])
    ids = ids[mask_random]
    ra = ra[mask_random]
    dec = dec[mask_random]
    z_true = z_true[mask_random]
    g1 = g1[mask_random]
    g2 = g2[mask_random]

# Load N(z)'s
case = 'test_mcal_bpzmof_unblind'
labels = [r'\textsc{Metacalibration}']
path_cats = '/Volumes/Data/y1_shear_tests/cats/'
path_runs = '/Volumes/Data/y1_shear_tests/runs/'
path_nofzs = path_runs+case+'/nofzs/'
nzs = np.loadtxt(path_nofzs + 'source')
nzs = nzs.T
nz_zmean = nzs[0]
nz_zlims = np.zeros(len(nz_zmean) + 1)
nz_zlims[:-1] = nz_zmean - (nz_zmean[2] - nz_zmean[1])/2. 
nz_zlims[-1] = nz_zmean[-1] + (nz_zmean[2] - nz_zmean[1])/2.
nz_all = []
for s in range(0,4):
    nz_all.append(nzs[s+1])

nz_all = np.array(nz_all)    
zbins, nz_weights = get_zbins_match_nz(zs = z_true, nz_zlims=nz_zlims, nzs=nz_all)


# Plot the distributions to make sure they match
colors = ["#B32222", '#ff9933', "#455A64", 'lightseagreen', 'saddlebrown', 'black']
fig, ax = plt.subplots()
for i in range(4):
    mask = zbins == i
    ax.hist(z_true[mask], bins=25, weights=nz_weights[mask], label = 'Bin %d'%(i+1), histtype='step', normed = True, lw = 2, color = colors[i])
    ax.plot(nz_zmean, nzs[i+1], color = colors[i], lw = 2)

ax.set_xlabel('Redshift')
ax.set_xlim(0, 1.6)
ax.legend()
plt.savefig('/Users/Judit/Dropbox/IFAE/y1/y1_shear_tests/histogram_mice_matched.pdf')


for i in range(4):
    print 'Saving file %d...'%i
    mask = zbins == i
    c1 = pf.Column(name='RA', format='E', array=ra[mask])
    c2 = pf.Column(name='DEC', format='E', array=dec[mask])
    c3 = pf.Column(name='ZTRUE', format='E', array=z_true[mask])
    c4 = pf.Column(name='W', format='E', array=nz_weights[mask])
    c5 = pf.Column(name='e1', format='E', array=g1[mask])
    c6 = pf.Column(name='e2', format='E', array=g2[mask])
    c7 = pf.Column(name='ID', format='E', array=ids[mask])

    CC = [c1,c2,c3,c4,c5,c6]
    hdu = pf.new_table(CC, nrows=len(ra[mask]))
    hdu.writeto('/Volumes/Data/y1_shear_tests/sims/mice/mice2_shear_fullsample_bin%d.fits'%(i+1), clobber=True)
