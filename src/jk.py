import numpy as np
import os
import kmeans_radec
#import pyfits as pf
import astropy.io.fits as pf
import matplotlib.pyplot as plt
import matplotlib.cm as cm
cmap = cm.Spectral
#from skymapper import *
from info import config

def jk_centers(ra_r, dec_r, path):
	np.random.seed(0)
	rand = np.random.random(len(ra_r))  
	maskr = rand>0.999
	ra_r_down = ra_r[maskr]
	dec_r_down = dec_r[maskr]

	#Filling radec array...
	radec_r = np.zeros((len(ra_r_down),2))
	radec_r[:,0] = ra_r_down
	radec_r[:,1] = dec_r_down

	#Starting kmeans...
	km = kmeans_radec.kmeans_sample(radec_r,config['njk'],maxiter=500,tol=1e-05) 

	#Did we converge?
	assert km.converged > 0, 'Kmeans did not converge! Try more iterations.'

	#Saving the kmeans centers
	np.savetxt('%s/jk_centers'%path,km.centers)

def jk(ra, dec, path):
	#Loading the jk centers
	fname = '%s/jk_centers'%path
	if os.path.isfile(fname): centers = np.loadtxt(fname)
	else: return 0

	#Filling radec array...
	radec = np.zeros((len(ra),2))
	radec[:,0] = ra
	radec[:,1] = dec

        # split in chuncks
        n = 10
	#Filling radec array...
        step = len(ra)/n

        jk = []
        for i in range(n+1):
                jki = kmeans_radec.find_nearest(radec[i*step:(i+1)*step], centers)
                jk.extend(jki)

        jk = np.array(jk)
        
        '''
	#Filling radec array...
	radec = np.zeros((len(ra),2))
	radec[:,0] = ra
	radec[:,1] = dec
	jk = kmeans_radec.find_nearest(radec, centers)
        '''        
	return jk


def plot_jk(paths, config, ra, dec, jk):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig = plt.figure(figsize=(6.,4.5))
	ax = fig.add_subplot(111, aspect='equal')
	proj = createConicMap(ax, ra, dec, proj_class=AlbersEqualAreaProjection, bgcolor='w')
	x,y = proj(ra, dec)
	marker = 's'
	markersize = getMarkerSizeToFill(fig, ax, x, y)
	sc = ax.scatter(x,y, c=jk, edgecolors='None', marker=marker, s=markersize, cmap=cmap, rasterized=True, zorder=1)
	ax.set_xlabel('RA',size='large')
	ax.set_ylabel('Dec',size='large')
	plt.savefig(paths['plots_config'] + 'jk.pdf')
	#plt.savefig(os.path.join(paths['redmagic',config['redmagic_v']],'jk.pdf'))

