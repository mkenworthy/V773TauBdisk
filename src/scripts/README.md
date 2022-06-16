# Scripts

Astrometry of the SPHERE data is carried out on `V773_SPHERE_median.fits` and uses the Jupyter notebook `measure_ABC_SPHERE_astrometry_with_emcee.ipynb`
.

Here we iterated manually a starting position for A and C in the cropped image.

Used simplex to fit for B, then did `emcee` to refine the values and locate the position of B and C relative to A.

Chains are stored in `astrom_ABC_chain.h5`

This chain is then used in `make_astrometry_ABC_plot.py` to make the figure `astrometry_ABC.pdf`