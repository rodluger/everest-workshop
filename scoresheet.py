import tools
import numpy as np
import matplotlib.pyplot as pl
from scipy.stats import norm
from truths import *

# Show answers?
show_answers = False

# Group 1 results
group1 = {
          1: [(0.0003, 0.0001)],
          2: [()],
          3: [()],
          4: [()],
          5: [()],
          6: [()],
          7: [()],
          8: [()],
          9: [()],
         10: [()]
         }

# Group 2 results
group2 = {
          1: [(0.0002, 0.0001)],
          2: [()],
          3: [()],
          4: [()],
          5: [()],
          6: [()],
          7: [()],
          8: [()],
          9: [()],
         10: [()]
         }

# Set up the figure
fig, ax = pl.subplots(2, 5, figsize = (12, 5))
ax = ax.flatten()
ax[1].set_title('De-trend THEN search', color = 'r', y = 1.1)
ax[2].set_title('vs.', color = 'k', y = 1.1)
ax[3].set_title('De-trend AND search', color = 'b', y = 1.1) 

for ID in range(1,11):

  # Plot group 1 results
  try:
    for depth1, std1 in group1[ID]:
      x = np.linspace(depth1 - 5 * std1, depth1 + 5 * std1)
      ax[ID-1].fill_between(x, np.zeros_like(x), norm.pdf(x, depth1, std1), color = 'r', alpha = 0.2)
  except ValueError:
    pass
  
  # Plot group2 results
  try:
    for depth2, std2 in group2[ID]:
      x = np.linspace(depth2 - 5 * std2, depth2 + 5 * std2)
      ax[ID-1].fill_between(x, np.zeros_like(x), norm.pdf(x, depth2, std2), color = 'b', alpha = 0.2)
  except ValueError:
    pass
  
  # Show the answers?
  if show_answers:  
    if t0[ID-1] is None:
      ax[ID-1].annotate('No planet!', xy = (0.5, 0.5), ha = 'center', va = 'center', xycoords = 'axes fraction', color = 'r', zorder = 99)
      ax[ID-1].set_xticks([])
      ax[ID-1].set_yticks([])
    else:
      if inject[ID - 1]:
        ax[ID-1].axvline(depth[ID - 1], color = 'k', alpha = 0.75)
  
  # Tweak appearance
  for tick in ax[ID-1].get_xticklabels():
    tick.set_fontsize(7)
  ax[ID-1].set_yticks([])
  ax[ID-1].annotate(ID, xy = (0.05, 0.95), xycoords = 'axes fraction', ha = 'center', va = 'center', fontweight = 'bold')


# Show!
pl.show()