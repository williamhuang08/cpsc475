import numpy as np
def gaborfilter(theta, wavelength, phase, sigma, aspect, ksize=None):

	"""
	GB = GABORFILTER(THETA, WAVELENGTH, PHASE, SIGMA, ASPECT, KSIZE)
	creates a Gabor filter GB with orientation THETA (in radians),
	wavelength WAVELENGTH (in pixels), phase offset PHASE (in radians),
	envelope standard deviation SIGMA, aspect ratio ASPECT, and dimensions
	KSIZE x KSIZE. KSIZE is an optional parameter, and if omitted default
	dimensions are selected.
	 """

	if ksize is None:
	 	ksize = 8*sigma*aspect


	if type(ksize) == int or len(ksize) == 1:
	 	ksize = [ksize, ksize]


	xmax = np.floor(ksize[1]/2.)
	xmin = -xmax
	ymax = np.floor(ksize[0]/2.)
	ymin = -ymax

	xs, ys = np.meshgrid(np.arange(xmin,xmax+1), np.arange(ymax,ymin-1,-1))

	# Your code here
