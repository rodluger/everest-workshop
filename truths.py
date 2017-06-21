#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
truths.py
---------

'''

# Target IDs
target = [211987102,  # 01
          211492412,  # 02
          211541097,  # 03
          211562654,  # 04
          211537731,  # 05
          211906824,  # 06
          212082447,  # 07
          211987102,  # 08
          211492412,  # 09
          211541097,  # 10
          ]

inject = [True,
          True,
          True,
          False,
          False,
          True,
          True,
          False,
          True,
          False,
         ]

# Time of first transit          
t0 =     [2310.873, 
          2307.631, 
          2308.812, 
          2314.782,   # Real planet (Pope et al. 2016)
          None,       # Not a planet host!
          2309.444, 
          2308.121, 
          None,       # Not a planet host!
          2315.86,
          None,       # Not a planet host!
          ]

# Period (days)  
per =    [15.789, 
          11.232, 
          10.1932, 
          10.792,     # Real planet (Pope et al. 2016)
          None,       # Not a planet host!
          14.723,
          13.383,
          None,       # Not a planet host!
          5.6375,
          None,       # Not a planet host!
          ]

# Duration (days)       
dur =    [0.2, 
          0.2, 
          0.2, 
          0.2,        # Real planet
          None,       # Not a planet host!
          0.2,
          0.2,
          None,       # Not a planet host!
          0.2,
          None,       # Not a planet host!
          ]

# Depth       
depth =  [0.0002, 
          0.001, 
          0.00015, 
          0.00068,     # Real planet (Pope et al. 2016)
          None,       # Not a planet host!
          0.0002,
          0.0003,
          None,       # Not a planet host!
          0.0004,
          None,       # Not a planet host!
          ]