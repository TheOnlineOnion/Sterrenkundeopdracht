import numpy as np
import astropy as ap
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.utils.data import download_file

MIN_HELDERHEID = 700

def verwerkte_afbeelding(path):
    """Voert een verwerkte afbeelding als np-array van pixelwaarden
    uit bij invoer van een pad naar een .fits-bestand.
    """

    # Pixelwaardes ophalen uit hdu_list van .fits-bestand
    hdu_list = fits.open(path) 
    image_data_raw = hdu_list[0].data
    hdu_list.close()

    # Afbeelding croppen, het is 0<=y<3330 en 0<=x<5363
    image_data = image_data_raw[0:3330, 0:5363]

    # Data invoeren in np-array
    stack = 1
    final_image = np.zeros_like(image_data)
    for n in range(0, stack):
        final_image += image_data

    return final_image

def maak_contour(img):
    """Voert een afbeelding van rode contouren rondom de sterren uit
    bij invoer van een verwerkte afbeelding (np-array).
    """
    cs = plt.contourf(img, levels=[MIN_HELDERHEID, MIN_HELDERHEID + 1], colors='red', extend='both', aspect='equal')
    cs.cmap.set_over('red')
    cs.cmap.set_under('blue')
    cs.changed()

    return cs

def maak_contourkaart(img):
    """Voert een bool array uit met 'True' waar in de contour rode
    pixels / sterren zijn.
    """
    bool_arr = img > MIN_HELDERHEID
    return bool_arr

def eet_rode_pixels(contour, x, y):  # x en y waardes van in te voeren pixel
    """Recursieve functie die rode pixels ('True') in blauwe pixels ('False') verandert,
    rondom de ingevoerde pixel.
    """
    if x < 0 or y < 0 or x >= contour.shape[0] or y >= contour.shape[1]:  # buiten randen
        return
    
    pixel = contour[x, y]
    if not pixel:  # pixel is niet 'rood'
        return
    contour[x, y] = False  # maak pixel 'blauw', zodat deze niet vaker wordt geteld.
    
    eet_rode_pixels(contour, x + 1, y)
    eet_rode_pixels(contour, x, y + 1)
    eet_rode_pixels(contour, x - 1, y)
    eet_rode_pixels(contour, x, y - 1)

    return

def tel_sterren(contour):
    """Voert het aantal groepen rode pixels, dus sterren,
    in de ingevoerde contourplot (als np-array) uit.
    """
    sterren = 0
    for x in range(contour.shape[0]):
        for y in range(contour.shape[1]):
            if contour[x, y]:
                eet_rode_pixels(contour, x, y)
                sterren += 1

    return sterren

def main(path):
    """Bij invoer van het pad naar een .fits-bestand: voert een plot uit
    van de afbeelding en bijbehorende contourplot, ofwel een histogram
    van helderheidswaardes.
    """

    modus = input('Voer "r" in voor het resultaat of "h" voor de histogram of "t" om sterren te tellen: ')

    afb = verwerkte_afbeelding(path)

    if modus == "r":

        # Maak een plot aan
        fig, axes = plt.subplots(2, 1)
        fig.set_figheight(6)
        fig.set_figwidth(6)

        contour = maak_contour(afb)

        # Verwerkte afbeelding en contourplot onder elkaar plotten
        # Merk op: vmin en vmax zijn de minimum en maximum lichtwaardes, zodat er geen ruis 
        # is in de afbeelding. Deze waarden zijn uit een histogram verkregen.
        axes[0].imshow(afb, cmap="grey", vmin=0.65e3, vmax=0.75e3, origin='lower', aspect='auto') 
        axes[1] = contour
    elif modus == "h":
        plt.hist(afb.flatten(), bins="auto")
    elif modus == "t":
        contourkaart = maak_contourkaart(afb)
        print(tel_sterren(contourkaart))

    plt.show()

    return

if __name__ == '__main__':
    # Voer in 'main' als argument het pad naar de .fit (of .fits) file in!!
    main(r"C:\Users\jelle\OneDrive - UvA\Dubba_wisnat\Sem 2 per 4\Sterrenkundeopdracht (AVT)\Light_Stack_19frames_20230906_TELESCOOP1_ELISE_5DOUBLECLUSTER_5sec_Bin1_filter-R_0.5C_gain0_2023-09-06_220632.fit")