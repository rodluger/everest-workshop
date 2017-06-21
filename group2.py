import tools
import numpy as np

# USER: The target number, the planet phase, and the planet period
ID = 1
phase = None
per = None

# Load the light curve
star = tools.Load(ID)

# Search and plot
if (phase is None) or (per is None):

  # Compute the delta-chi squared metric.
  # You can tune the periods and phases kwargs for a more 
  # constrained/higher resolution search grid.
  star.search(joint_fit = True, 
              periods = np.linspace(5., 20., 1000), 
              phases = np.linspace(0., 1., 1000))

else:

  # Once you know the phase and the period, compute the depth posterior
  depth, depth_uncert = star.compute_depth(joint_fit = True, phase = phase, per = per)

  # Once you know the depth, plot the folded light curve!
  star.plot_folded(phase = phase, per = per, depth = depth)

  # Print the result
  print("Depth = %.6f +/- %.6f" % (depth, depth_uncert))