import numpy as np

import astropy as ap
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats, SigmaClip

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.utils.data import download_file
from astropy.visualization import ImageNormalize, AsymmetricPercentileInterval

import photutils as pu
from photutils.aperture import CircularAperture
from photutils.background import Background2D, MeanBackground
from photutils.detection import DAOStarFinder

# Je moet hieronder zelf handmatig de file location veranderen naar waar jij je fits files hebt opgeslagen
path = "C:/Users/MeintPostSecuMailerC/Documents/vnv_natuurkunde/Natuurkunde Opdracht/FITS images/Light_Stack_19frames_20230906_TELESCOOP1_ELISE_5DOUBLECLUSTER_5sec_Bin1_filter-R_0.5C_gain0_2023-09-06_220632.fit"

# We maken er een CCD datatypen van, want daar staat de header en image data beide in
ccd_image = CCDData.read(path, unit="adu")
image_data = ccd_image.data
header = ccd_image.header

norm_image = ImageNormalize(ccd_image, interval=AsymmetricPercentileInterval(30, 99.5))

# 2-D background estimation
sigma_clip = SigmaClip(sigma=3., maxiters=5)
bkg_estimator = MeanBackground()
bkg = Background2D(ccd_image, box_size=200, filter_size=(11, 11), 
                   mask=ccd_image.mask, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

# Calculate the 2D background subtraction, maintaining metadata, unit, and mask
ccd_2d_bkgdsub = ccd_image.subtract(bkg.background)

# Set up the figure with subplots
fig, ax1 = plt.subplots(1, 1, figsize=(8, 8), sharey=True)
plt.tight_layout()


ax1.set_ylabel('Y (pixels)')
ax1.set_xlabel('X (pixels)')

# Plot the 2D-subtracted data
norm_image_sub = ImageNormalize(ccd_2d_bkgdsub, interval=AsymmetricPercentileInterval(30, 99.5))

fitsplot = ax1.imshow(np.ma.masked_where(ccd_2d_bkgdsub.mask, ccd_2d_bkgdsub), norm=norm_image_sub, cmap='viridis')
ax1.set_title('2D Background-Subtracted Data')

# Define the colorbar...
cbar = fig.colorbar(fitsplot, location='right', shrink=0.6)

cbar.set_label(r'Counts ({})'.format(ccd_image.unit.to_string('latex')),
               rotation=270, labelpad=30)

plt.show()