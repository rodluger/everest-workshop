#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
generate.py
----------

planethost01.fits -- AWESOME


'''

from __future__ import division, print_function, absolute_import, unicode_literals
import everest
import os
import logging
log = logging.getLogger(__name__)

# Injected planet params
target = [211987102, 211492412, 211541097, 211733267]
t0 = [2310.873, 2307.631, 2308.812, None]
per = [15.789, 11.232, 10.1932, None]
dur = [0.2, 0.2, 0.2, None]
depth = [0.0002, 0.001, 0.0002, None]

# Loop over all targets
for i in range(len(target)):
  
  log.info("Generating light curve for target %d/%d..." % (i + 1, len(target)))
  
  # De-trend the original light curve
  star_orig = everest.rPLD(target[i], debug = True)
  
  # Injection or real planet?
  if t0 is not None:
  
    # Inject the planet, de-trend, and generate a FITS file
    star = everest.Inject(target[i], inj_model = 'rPLD', t0 = t0[i], per = per[i], dur = dur[i], 
                          depth = depth[i], mask = False, debug = True)
    star.cdppg = 0
    star.clobber = True
    everest.fits.MakeFITS(star, fitsfile = 'tmp.fits')
    os.rename(os.path.join(star.dir, 'tmp.fits'), 'planethost%02d.fits' % (i + 1))
  
  else:
    
    # Real planet!
    everest.fits.MakeFITS(star_orig, fitsfile = 'tmp.fits')
    os.rename(os.path.join(star_orig.dir, 'tmp.fits'), 'planethost%02d.fits' % (i + 1))