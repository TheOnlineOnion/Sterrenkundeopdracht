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

# We gebruiken de sigma-clip functionaliteit om een 2 dimensionale waarde plot te maken van de achtergrond, zodat we deze er af kunnen halen en overal een locale gemiddelde achtergrond waarde van 0 kunnen bereiken
sigma_clip = SigmaClip(sigma=3., maxiters=5)
bkg_estimator = MeanBackground()
bkg = Background2D(ccd_image, box_size=200, filter_size=(11, 11), 
                   mask=ccd_image.mask, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

# We verwijderen nu de locale gemmidelde achtergrond waarde van elke pixel die niet een ster is, zodat de steren better gedetecteerd kunnen worden
ccd_2d_bkgdsubtr = ccd_image.subtract(bkg.background)

mean, median, std = sigma_clipped_stats(ccd_2d_bkgdsubtr.data, sigma=3, maxiters=5)

# We maaken een mask array aan
mask = np.zeros_like(ccd_2d_bkgdsubtr, dtype=bool)

# en we maskeren de buitenste 10 pixels, aangezien die geen echte data representeren op de plot
mask[:10, :] = mask[-10:, :] = mask[:, :10] = mask[:, -10:] = 1.0

# We maken een object aan van de DAOStarFinder classe, met een full width half maximum van 6.2
fwhm = 6.2
daofinder = DAOStarFinder(threshold=6 * std, fwhm=fwhm)

# Zoek nu de sterren, zonder de buitenlaag mee te nemen
sources = daofinder(ccd_2d_bkgdsubtr.data, mask=mask)
print(f"Found {len(sources)} sources in the image")

# We plaatsen cirkels rondom de x en y waardes van de sterren, om te zien welke sterren de starfinder wel en niet pakt
positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
apertures = CircularAperture(positions, r=3 * fwhm)

# De afbeelding plotten, met colorbar
fig, ax1 = plt.subplots(1, 1, figsize=(8, 8), sharey=True)
plt.tight_layout()

# Labels voor de assen aanmaken
ax1.set_ylabel('Y (pixels)')
ax1.set_xlabel('X (pixels)')

# De data zonder achtergrond nomralizeren en daarna plotten
norm_image_sub = ImageNormalize(ccd_2d_bkgdsubtr, interval=AsymmetricPercentileInterval(30, 99.5))

fitsplot = ax1.imshow(np.ma.masked_where(ccd_2d_bkgdsubtr.mask, ccd_2d_bkgdsubtr), norm=norm_image_sub, cmap='viridis')
ax1.set_title('2D Background-Subtracted Data')

# De cirkels plotten
apertures.plot(color="red")

# De colorbar aanmaken en op de plot plaatsen
cbar = fig.colorbar(fitsplot, location='right', shrink=0.6)

cbar.set_label(r'Counts ({})'.format(ccd_image.unit.to_string('latex')),
               rotation=270, labelpad=30)

plt.show()