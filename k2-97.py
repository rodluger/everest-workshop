from everest import Everest, TransitModel
from everest.transit import Get_RpRs
import matplotlib.pyplot as pl

# The true planet properties for K2-97 from
# Pope et al. (2016), hosted at the NASA Exoplanet Archive 
# (exoplanetarchive.ipac.caltech.edu)
EPIC = 211351816
per = 8.406
t0 = 2309.05
rhos = 0.4
RpRs = 0.023577

# Load the EVEREST light curve 
star = Everest(EPIC)

# Show the EVEREST data validation summary
star.dvs()

# Plot the light curve interactively
star.plot()

# Conduct a transit search. We specify the **approximate** planet period,
# which constrains the transit duration. In a full search, we'd have to
# loop over a (coarse) period grid.
time, depth, vardepth, delchisq = star.search(rhos = rhos, per = 10., clobber = False)

# Plot the results, and indicate the actual transits
fig, ax = pl.subplots(1, figsize = (10, 4))
ax.plot(time, delchisq, lw = 1)
ax.set_ylabel(r'$\Delta \chi^2$', fontsize = 18)
ax.set_xlabel('Time (days)', fontsize = 18)
ax.set_xlim(time[0], time[-1])
for n in range(9):
  ax.axvline(t0 + n * per, color = 'r', ls = '--', alpha = 0.5, zorder = -1)
pl.show()

# Compute the joint transit model and calculate the transit depth
star.transit_model = TransitModel('b', per = per, t0 = t0, rhos = rhos)
star.plot_transit_model(fold = 'b')

# Conver to a radius ratio
RpRs_everest = Get_RpRs(star.transit_depth, per = per, rhos = rhos)
print("EVEREST transit depth: %.3e" % RpRs_everest)