import tools
import numpy as np
  
# USER: The target number
ID = 1

# USER: Perform a joint instrumental/systematics fit?
joint_fit = True

# Load the light curve
star = tools.Load(ID)

# Compute the delta-chi squared metric: search on de-trended data
time, ml_depth, depth_variance, delta_chisq = star.search(joint_fit = joint_fit)