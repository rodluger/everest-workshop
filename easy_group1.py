import tools
import numpy as np
  
# USER: The target number
ID = 4

# Load the light curve
star = tools.Load(ID)

# Search!
star.search()