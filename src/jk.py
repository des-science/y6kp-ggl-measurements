import numpy as np
import os
import kmeans_radec
import pyfits as pf
import matplotlib.pyplot as plt
import matplotlib.cm as cm
cmap = cm.Spectral
from skymapper import *

def jk_centers(paths, config):
	#Loading random subsample for kmeans generation 
	data = pf.getdata(paths['redmagic',config['redmagic_v'],'hidens_randoms'])
	ra_r = data['ra']
	dec_r = data['dec']

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
	np.savetxt(os.path.join(paths['redmagic',config['redmagic_v']],'jk_centers'),km.centers)


def jk(paths, config, ra, dec):
	#Loading the jk centers
	fname = os.path.join(paths['redmagic',config['redmagic_v']],'jk_centers')
	if os.path.isfile(fname): centers = np.loadtxt(fname)
	else: jk_centers(paths,config)

	#Filling radec array...
	radec = np.zeros((len(ra),2))
	radec[:,0] = ra
	radec[:,1] = dec

	jk = kmeans_radec.find_nearest(radec, centers)
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
