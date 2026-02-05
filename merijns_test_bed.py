#%% Afbeelding
import numpy as np
import astropy as ap
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.utils.data import download_file

# Je moet hieronder zelf handmatig de file location veranderen naar waar jij je fits files hebt opgeslagen
path = "C:/Users/MeintPostSecuMailerC/Documents/vnv_natuurkunde/Natuurkunde Opdracht/FITS images/Dubbelcluster_2.fit"

informatie = fits.getheader(path)
print(informatie)

hdu_list = fits.open(path) 
image_data_l = hdu_list[0].data
hdu_list.close()

# Afbeelding croppen, het is 0<=y<3330 en 0<=x<05363
image_data = image_data_l[0:3330, 0:5363]
image_data = image_data*0.75

# l is hoe vaak je de afbeelding wil stacken
l = 1
final_image = np.zeros_like(image_data)
for n in range(0, l):
    final_image += image_data
print(f"This is stacked {l} times")

#Histogram om vmin en vmax te verkrijgen
#image_hist = plt.hist(final_image.flatten(), bins="auto")

fig = plt.figure()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#vmin en vmax verkegen uit de histogram
ax.imshow(final_image, cmap="grey", vmin=1.5e3, vmax=2.1e3) 
#ax.contour(final_image, levels=10, colors='red', linewidths=1.5)

plt.show()

