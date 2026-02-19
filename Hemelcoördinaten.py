from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

hdu = fits.open("new-image.fits")[0]

w = WCS(hdu.header)

ny, nx = hdu.data.shape

y, x = np.mgrid[0:ny, 0:nx]

ra, dec = w.all_pix2world(x, y, 0)

print("RA shape:", ra.shape)
print("DEC shape:", dec.shape)
print(f'de pixel (100, 100) heeft hemelcoördinaten: ra = {ra[100][100]} & dec = {dec[100][100]}')