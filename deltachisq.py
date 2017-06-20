import tools
import numpy as np
  
# USER: The target number, and whether or not to do a joint instrumental/systematics fit
ID = 3
joint_fit = False

# Load the light curve
star = tools.Load(ID)

# Compute the delta-chi squared metric: search on de-trended data
time, ml_depth, depth_variance, delta_chisq = star.search(joint_fit = joint_fit)