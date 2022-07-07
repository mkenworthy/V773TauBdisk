import paths
import numpy as np
import matplotlib.pyplot as plt
#from rotscatx import *
from astropy.io import fits
import corner
import emcee

# this uses chains generated by measure_ABC_SPHERE_astrometry_with_emcee.ipynb

def sub_stars(Ax, Ay, dx, dy, df, im):
    Arot = rotscatx(im, (-Ay,-Ax), 180., 1., ([0,0]), 2)
    # now move C to the point of B...
    B_pos = rotscatx(im, (0,0), 0., 1., ([(-dy),(-dx)]), 2)*df

    return (im-Arot-B_pos)

from scipy import ndimage

def rotscatx(im, p0, rot=10., sca=1., p1=([0.,0.]), order = 1 ):
    """rotate, scale and translate using one affine transformation

    rotate and scale a 2D image with optional offset

    :Parameters:
      - `im` (float) - a 2D image
      - `p0` (float) - centre point for rotation and scaling
      - `rot` (float) - rotation angle (degrees)
      - `sca` (float) - scaling
      - `p1` (float) - translation vector after rotation/scaling

    :Returns:
      - `im2` (float) - the transformed image

    :Examples:

    rotate img by 30 degrees around point 10,10 and translate by 2,2
    interpolate using splines of order 2
    >>> rotscatx(img, (-10,-10), 30., 1., ([2.,2.]), 2)

    """

    # affine discussion based on:
    # http://answerpot.com/showthread.php?3899770-%231736%3A+Proper+order+of+translation+and+rotation+in+scipy.ndimage.affine_transform

    # make a rotation matrix
    theta = np.deg2rad(rot)
    rotate = np.identity(2)

    rotate[:2,:2] = [np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]

    # make a scaling matrix
    scale = np.identity(2)
    scale *= sca

    # make affine matrix
    affine = np.dot(scale, rotate)

    # affine transform the input centre point
    p0aff = np.dot(affine, p0)

    return ndimage.interpolation.affine_transform(im, affine, offset=p0aff-p0+p1, order=order)


# read in SPHERE image
image_file = 'V773_SPHERE_median.fits'
image_data = fits.getdata(paths.data / image_file, ext=0)

# crop out star from the image
xc,yc = (512,512)
b = 55
imc = image_data[yc-b:yc+b,xc-b:xc+b]

# positions of A and C in the cropped image
C_y, C_x = (38.087, 40.694)
x_A, y_A = (54.495,54.489)

# chain of walkers from astrometric fitting with emcee
file_chain = paths.data / 'astrom_ABC_chain.h5'

reader = emcee.backends.HDFBackend(file_chain,read_only=True)

ndim, nwalkers = 3, 100
labels = ["$dx$","$dy$", "$f$"]

chain = reader.get_chain()
lnprob = reader.get_log_prob()

# see what the chains look like, skip a burn in period if desired
burn = 1000

# samples for deriving fitted parameters and uncertainties

#samples = chain[burn:,].reshape((-1, ndim))
samples = chain[burn:,].reshape((-1, ndim))

fig, ax = plt.subplots(ndim,ndim, figsize=(8,8))

fig2 = corner.corner(samples, labels=labels, show_titles=False, fig=fig)
fig.savefig("_check_astrometry_ABC_corner_plot.pdf")

for i in range(ndim):
        mcmc = np.percentile(samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
        txt = txt.format(mcmc[1], q[0], q[1], labels[i])
        print(txt)

dx = np.percentile(samples[:, 0], [16, 50, 84])[1]
dy = np.percentile(samples[:, 1], [16, 50, 84])[1]
f = np.percentile(samples[:, 2], [16, 50, 84])[1]

fig, (ax2, ax3) = plt.subplots(1,2,figsize=(6,3))

Arot = rotscatx(imc, (-y_A,-x_A), 180., 1., ([0,0]), 2)

ax2.imshow(imc-Arot,
           origin='lower',
           vmin=-100,vmax=100
          )

model = sub_stars(x_A, y_A, dx, dy, f, imc)

ax3.imshow(model,
           origin='lower',
           vmin=-100,vmax=100
          )

ax2.scatter(x_A, y_A,color='red')
ax3.scatter(C_x+15.9, C_y+22.05,color='red',alpha=0.5)

fig.savefig(paths.figures / "astrometry_ABC.pdf", bbox_inches='tight')
