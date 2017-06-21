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

# Once you know the phase and the period, compute the depth posterior
# depth, depth_uncert = star.compute_depth(phase = 0., per = 10.)

# Once you know the depth, plot the folded light curve!
# star.plot_folded(phase = 0., per = 10., depth = depth)
