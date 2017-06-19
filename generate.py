#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
generate.py
----------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import everest
import os
import logging
from truths import *
log = logging.getLogger(__name__)

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
    os.rename(os.path.join(star.dir, 'tmp.fits'), 'target%02d.fits' % (i + 1))
  
  else:
    
    # Real planet!
    everest.fits.MakeFITS(star_orig, fitsfile = 'tmp.fits')
    os.rename(os.path.join(star_orig.dir, 'tmp.fits'), 'target%02d.fits' % (i + 1))