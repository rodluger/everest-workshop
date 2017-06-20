import tools
import numpy as np
  
# USER: The target number
ID = 1

# Load the light curve
star = tools.Load(ID)

# Compute the delta-chi squared metric.
# This function returns a time array, the maximum likelihood depth
# array, the variance on the depth estimate, and the delta chi squared
# array, in case you want to play around with them. 
# You can also tune the periods and phases kwargs for a more 
# constrained/higher resolution search grid.
time, ml_depth, depth_variance, delta_chisq =\
  star.search(joint_fit = True, 
              periods = np.linspace(5., 20., 1000), 
              phases = np.linspace(0., 1., 1000))


# Once you figured out the phase, period, and depth of your planet,
# you can plot the folded light curve
phase = 0.2695
per = 15.789
depth = 0.0002
t0 = star.time[0] + phase * per

# You can compute the light curve de-trended with
# a joint instrumental/transit model fit and plot the folded transit
star.transit_model = tools.TransitModel(star.time, t0 = t0, per = per, depth = depth)
star.compute_joint()
depth = star.transit_depth

# Print to screen and plot!
print("PERIOD:        %.3f days" % per)
print("FIRST TRANSIT: %.3f days" % t0)
print("DEPTH:         %.3e" % depth)
star.plot_folded(t0, per, depth)