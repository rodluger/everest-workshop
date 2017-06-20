import tools
import numpy as np

# USER: The target number
ID = 1

# Load the light curve
star = tools.Load(ID)

# Compute the delta-chi squared metric.
# You can tune the periods and phases kwargs for a more 
# constrained/higher resolution search grid.
star.search(joint_fit = True, 
            periods = np.linspace(5., 20., 1000), 
            phases = np.linspace(0., 1., 1000))

# Once you know the phase, period, and depth,
# plot the folded light curve!
#
# star.plot(joint_fit = True, phase = 0., per = 10., depth = 0.001)