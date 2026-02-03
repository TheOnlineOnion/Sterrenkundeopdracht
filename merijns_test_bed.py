import numpy as np
import astropy as ap
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.utils.data import download_file

hdu_list = fits.open("C:/Users/MeintPostSecuMailerC/Documents/vnv_natuurkunde/Natuurkunde Opdracht/FITS images/Light_Stack_275frames_dubbelcluster_5sec_Bin2_filter-L_1.0C_gain111_2025-01-13_221943.fit") # Je moet hier zelf handmatig de file location veranderen naar waar jij je fits files hebt opgeslagen
hdu_list.info()
image_data = hdu_list[0].data
hdu_list.close()
l = 1
final_image = np.zeros_like(image_data)
for n in range(0, l):
    final_image += image_data
print(f"This is stacked {l} times")

#image_hist = plt.hist(final_image.flatten(), bins="auto")
plt.imshow(final_image, cmap="grey", vmin=1.8e3, vmax=2.3e3) #vmin en vmax verkegen uit de histogram
plt.colorbar()

plt.show()