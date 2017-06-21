import tools
import numpy as np
  
# USER: The target number
ID = 1

# Load the light curve
star = tools.Load(ID)

# Compute the delta-chi squared metric.
# You can tune the periods and phases kwargs for a more 
# constrained/higher resolution search grid.
star.search(periods = np.linspace(5., 20., 1000), 
            phases = np.linspace(0., 1., 1000))

# Once you know the phase and the period, enter them here:
phase = None
per = None

if not (phase is None) and not (per is None):

  # Once you know the phase and the period, compute the depth posterior
  depth, depth_uncert = star.compute_depth(phase = phase, per = per)

  # Once you know the depth, plot the folded light curve!
  star.plot_folded(phase = phase, per = per, depth = depth)

  # Print the result
  print("Depth = %.6f +/- %.6f" % (depth, depth_uncert))