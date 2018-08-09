import numpy as np
from astropy.io import fits
from astrodendro import Dendrogram, pp_catalog

# Directory where continuum image is located
directory = '../../Data/Continuum/12m_7m_SD/'
# Continuum image
image = directory+'brick_1mm_continuum.fits'
# Output directory for data products
outputdir = '../Output_data_products/Dendro/'
outputdendro = outputdir+'brick_1mm_continuum_dendro.fits'

# read in image
hdu = fits.open(image)
header = hdu[0].header
data = hdu[0].data
hdu.close()
data = np.squeeze(data)
print np.shape(data)

# Image information
cont_rms        = 5e-4 # 0.5mJy/beam
beamsize_major  = header['BMAJ'] * 60.0 * 60.0 # arcsec
beamsize_minor  = header['BMIN'] * 60.0 * 60.0 # arcsec
beam            = [beamsize_major/206205., beamsize_minor/206205.]
beamsize        = np.sqrt([beamsize_major*beamsize_minor]) # arcsec
pixel_size      = np.abs(header['CDELT1'] * 60.0 * 60.0) # arcsec
beamarea        = beamsize_major * beamsize_minor * ((2.*np.pi)/(8.*np.log(2.)))
pixelarea       = pixel_size**2.0
npix            = beamarea / pixelarea

# Compute the dendrogram

min_value = 3.0*cont_rms
min_delta = 3.0*cont_rms
min_npix  = npix

# Toggle for computation
dendro    = Dendrogram.compute(data, min_value=min_value, min_delta=min_delta, min_npix=min_npix, verbose=True)
dendro.save_to(outputdendro)

dendro.viewer() # Toggle to view the dendrogram
