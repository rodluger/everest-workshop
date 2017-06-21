import tools
import numpy as np
import matplotlib.pyplot as pl
from scipy.stats import norm
from truths import *

fig, ax = pl.subplots(2, 5, figsize = (12, 5))
ax = ax.flatten()

for ID in range(1,11):
  if t0[ID - 1] is not None:
    star = tools.Load(ID)
    phase = (t0[ID - 1] - star.time[0]) / per[ID - 1]
    depth1, std1 = star.compute_depth(phase = phase, per = per[ID - 1], joint_fit = False)
    depth2, std2 = star.compute_depth(phase = phase, per = per[ID - 1], joint_fit = True)
    x = np.linspace(depth1 - 5 * std1, depth1 + 5 * std1)
    ax[ID-1].fill_between(x, np.zeros_like(x), norm.pdf(x, depth1, std1), color = 'r', alpha = 0.3)
    x = np.linspace(depth2 - 5 * std2, depth2 + 5 * std2)
    ax[ID-1].fill_between(x, np.zeros_like(x), norm.pdf(x, depth2, std2), color = 'b', alpha = 0.3)
    if inject[ID - 1]:
      ax[ID-1].axvline(depth[ID - 1], color = 'k', alpha = 0.75)
    for tick in ax[ID-1].get_xticklabels():
      tick.set_fontsize(7)
    ax[ID-1].set_yticks([])
  else:
    ax[ID-1].annotate('No planet!', xy = (0.5, 0.5), ha = 'center', va = 'center', xycoords = 'axes fraction', color = 'r')
    ax[ID-1].set_xticks([])
    ax[ID-1].set_yticks([])
  ax[ID-1].annotate(ID, xy = (0.05, 0.95), xycoords = 'axes fraction', ha = 'center', va = 'center', fontweight = 'bold')

ax[1].set_title('De-trend THEN search', color = 'r', y = 1.1)
ax[2].set_title('vs.', color = 'k', y = 1.1)
ax[3].set_title('De-trend AND search', color = 'b', y = 1.1) 

pl.show()