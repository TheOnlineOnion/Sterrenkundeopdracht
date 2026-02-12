import numpy as np
import matplotlib.pyplot as plt
import astropy as ap
import astropy.io.fits as fits
import os 
from astropy.visualization import astropy_mpl_style

intensity_increase = 10
increase_boundary = 2000

hdul = fits.open('C:/Users/bosel/OneDrive/Documenten/Light_Stack_275frames_dubbelcluster_5sec_Bin2_filter-L_1.0C_gain111_2025-01-13_221943 (2).fit')

#if 'hh.fits' not in os.listdir():
    #hdul.writeto('hh.fits')

plt.style.use(astropy_mpl_style)

image_data = fits.getdata('hh.fits')

image_data = np.where(image_data > increase_boundary, intensity_increase * image_data, image_data)

plt.imshow(image_data, cmap='gray')
plt.colorbar()
plt.show()