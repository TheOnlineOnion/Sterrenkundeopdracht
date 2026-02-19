import numpy as np

import astropy as ap
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats, SigmaClip

import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, AsymmetricPercentileInterval
import matplotlib.tri as tri

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
fwhm = 6.7
daofinder = DAOStarFinder(threshold=5.8 * std, fwhm=fwhm)

# Zoek nu de sterren, zonder de buitenlaag mee te nemen
sources = daofinder(ccd_2d_bkgdsubtr.data, mask=mask)
print(f"Found {len(sources)} sources in the image")

# We plaatsen cirkels rondom de x en y waardes van de sterren, om te zien welke sterren de starfinder wel en niet pakt
positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
apertures = CircularAperture(positions, r=3 * fwhm)


# De bevolkingsdichtheids diagram maken. We hebben een afbeelding van dy=3672 en dx = 5496. We maken blokjes van 102 hoog en 102 breed zodat we 72 blokjes verticaal en 107 blokjes horizontaal hebben
size_x = 102
size_y = 102

len_x = 5496
len_y = 3672

image_density_map = np.zeros_like(ccd_2d_bkgdsubtr)
small_image_density_map = np.zeros([int(len_y/size_y), int(len_x/size_x)])
for y in range(1, int(len_y/size_y)+1):
    for x in range(1, int(len_x/size_x)+1):
        num_stars = 0
        for stars in positions:
            if (x - 1) * size_x < (stars[0]) <= x * size_x and (y - 1) * size_y < (stars[1]) <= y * size_y:
                num_stars += 1
        image_density_map[(y - 1) * size_y:y * size_y, (x - 1) * size_x:x * size_x] = num_stars
        small_image_density_map[(y-1):y, (x-1):x] = num_stars

x_stars = np.linspace(0, int(len_x/size_x), int(len_x/size_x)+1)
y_stars = np.linspace(0, int(len_y/size_y), int(len_y/size_y)+1)
xv_stars, yv_stars = np.meshgrid(x_stars, y_stars)
x_stars_flat = np.hstack(xv_stars)
y_stars_flat = np.hstack(yv_stars)
z_stars = small_image_density_map.flatten()

stars_non_zero = np.nonzero(small_image_density_map)

x_final = np.array(stars_non_zero[0])
y_final = np.array(stars_non_zero[1])

z_final = [z_stars[x] for x in np.nonzero(z_stars)][0]

print(x_final, y_final, z_final)


# De afbeelding plotten, met colorbar
fig, axes = plt.subplots(1, 3, figsize=(14, 3))
fig.tight_layout(w_pad=5, pad=1)

# Labels voor de assen aanmaken
axes[0].set_ylabel('Y (pixels)')
axes[0].set_xlabel('X (pixels)')

# De data zonder achtergrond nomralizeren en daarna plotten
norm_image_sub = ImageNormalize(ccd_2d_bkgdsubtr, interval=AsymmetricPercentileInterval(30, 99.5))

fitsplot = axes[0].imshow(np.ma.masked_where(ccd_2d_bkgdsubtr.mask, ccd_2d_bkgdsubtr), norm=norm_image_sub, cmap='viridis')
axes[0].set_title('Achtergrond Verwijderde Data, met ster detectie cirkels')

# De cirkels plotten
apertures.plot(axes[0], color="red")

# De colorbar aanmaken en op de plot plaatsen
cbar = fig.colorbar(fitsplot, location='right', shrink=0.6)

cbar.set_label(r'Counts ({})'.format(ccd_image.unit.to_string('latex')),
               rotation=270, labelpad=30)

# De bevolkingsdichtheids diagram plotten
density = axes[1].imshow(image_density_map, cmap='viridis')
axes[1].set_title('Bevolkingsdichtheid kaart dubbelclustert')
axes[1].set_ylabel('Y (pixels)')
axes[1].set_xlabel('X (pixels)')

# colorbar tweede figuur
cbar2 = fig.colorbar(density, location='right', shrink=0.6)
cbar2.set_label('# Stars')

#axes[2].plot( y_final, x_final, 'ko', ms=3)
contour = axes[2].tricontourf(y_final, x_final, z_final, levels=[0, 1.75, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], cmap='RdBu_r', antialiased=False)
axes[2].set(xlim=(0, int(len_x/size_x)-1), ylim=(int(len_y/size_y)-1, 0))
axes[2].set_title('Bevolkingsdichtheids contour plot dubbelcuster')
axes[2].set_ylabel('Y (pixels)')
axes[2].set_xlabel('X (pixels)')

# colorbar tweede figuur
cbar3 = fig.colorbar(contour, location='right', shrink=0.6)
cbar3.set_label('# Stars')

plt.show()