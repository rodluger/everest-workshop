import everest
import numpy as np
import matplotlib.pyplot as pl

# Load everest
star = everest.Everest(201367065)

# Let's look at what we're dealing with
# There are three transiting planets buried in here.
pl.plot(star.time, star.fraw, 'k.', alpha = 0.3, ms = 2)
pl.show()

# Let's remove the really bad outliers
cut = np.where(star.fraw < 355000)
time = np.delete(star.time, cut)
fpix = np.delete(star.fpix, cut, axis = 0)
ntime, npix = fpix.shape

# Let's look at what we're dealing with
for n in range(npix):
  pl.plot(time, fpix[:,n])
 pl.show()

# Construct our design matrix
total_flux = np.sum(fpix, axis = 1).reshape(-1, 1)

# Take 1: First order PLD design matrix 
X = fpix / total_flux.reshape(-1, 1)

# Let's look at what we're dealing with
for n in range(npix):
 pl.plot(time, X[:,n])
 pl.show()

# De-trend!
w = np.linalg.solve(np.dot(X.T, X), np.dot(X.T, total_flux))
model = np.dot(X, w)
detrended_flux = total_flux - model

# Normalize it
detrended_flux += np.nanmedian(total_flux)

# Plot!
fig, ax = pl.subplots(2)
ax[0].plot(time, total_flux, 'k.', alpha = 0.3, ms = 2)
ax[0].plot(time, model, 'r-', alpha = 0.5)
ax[1].plot(time, detrended_flux, 'k.', alpha = 0.3, ms = 2)
pl.show()

# - # - # - # - # - # - # - #

# Take 2: Add a polynomial baseline
X = np.hstack((X, np.array([np.linspace(0, 1, ntime) ** n for n in range(50)]).T))

# De-trend!
w = np.linalg.solve(np.dot(X.T, X), np.dot(X.T, total_flux))
model = np.dot(X, w)
detrended_flux = total_flux - model

# Normalize it
detrended_flux += np.nanmedian(total_flux)

# Plot!
fig, ax = pl.subplots(2)
ax[0].plot(time, total_flux, 'k.', alpha = 0.3, ms = 2)
ax[0].plot(time, model, 'r-')
ax[1].plot(time, detrended_flux, 'k.', alpha = 0.3, ms = 2)
pl.show()