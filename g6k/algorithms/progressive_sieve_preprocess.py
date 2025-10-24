# -*- coding: utf-8 -*-
####
#
#   Copyright (C) 2018-2021 Team G6K
#
#   This file is part of G6K. G6K is free software:
#   you can redistribute it and/or modify it under the terms of the
#   GNU General Public License as published by the Free Software Foundation,
#   either version 2 of the License, or (at your option) any later version.
#
#   G6K is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with G6K. If not, see <http://www.gnu.org/licenses/>.
#
####



from __future__ import absolute_import
from __future__ import print_function
import sys
import time
from math import log
from g6k.siever import SaturationError
import logging



def norm(vec):
    return sum([_**2 for _ in vec])

def print_pro_randslicer_state(progressive_sieve):
    progressive_sieve.minl = min(progressive_sieve.g6k.l, progressive_sieve.minl)
    if progressive_sieve.phase != "down":
        print("\r %3d: ↑%3d      " % (progressive_sieve.g6k.r-progressive_sieve.l, progressive_sieve.g6k.r-progressive_sieve.g6k.l), end=' ')
    else:
        print("\r %3d: ↑%3d ↓%3d " % (progressive_sieve.g6k.r-progressive_sieve.l, progressive_sieve.g6k.r-progressive_sieve.minl, progressive_sieve.g6k.r-progressive_sieve.g6k.l), end=' ')
    sys.stdout.flush()


# Sieve (switching to gauss if needed) and test loop-breaking conditions
# Return true if we should continue
def wrapped_sieve(progressive_sieve):
    if progressive_sieve.phase == "init":
        alg = "gauss"
    else:
        alg = None

    cont = True
    try:
        with progressive_sieve.g6k.temp_params(saturation_ratio=progressive_sieve.g6k.params.saturation_ratio * progressive_sieve.sat_factor):
           
            # Match lifting effort to insertion strategy
            progressive_sieve.g6k(alg=alg)
            
            
    except SaturationError as e:
        if progressive_sieve.saturation_error == "skip":
            progressive_sieve.down_sieve = False
            logging.info("saturation issue: breaking progressive_sieve.")
            cont = False
        elif progressive_sieve.saturation_error == "weaken":
            logging.info("saturation issue: weakening progressive_sieve.")
            progressive_sieve.sat_factor /= 2.
        elif progressive_sieve.saturation_error == "ignore":
            pass
        else:
            raise e
    return cont




def progressive_sieve(g6k, kappa, blocksize,  dim4free,      # Main parameters
         start_up_n=50,     verbose=False,  saturation_error="weaken",                                                    
         ):
    """
    Run the progressive_sieve algorithm.

    :param g6k: The g6k object to work with
    :param c: target vector
    :param kappa: beginning of the block
    :param blocksize: dimension of the block (r=kappa+blocksize)
    :param dim4free: number of ``dimension for free'' [Ducas, Eurcrypt 2018]: Sieve context [l,r] where l=kappa+dim4free
    :param verbose: print progressive_sieve steps on the standard output.
    :param saturation_error: determines the behavior of pump when encountering saturation issue {"weaken",
        "skip", "ignore", "forward"}
    
    Generate sieve list for slicer preprocess.

    """
    progressive_sieve.l = kappa+dim4free  # noqa
    progressive_sieve.r = kappa+blocksize
    g6k.shrink_db(0)
    g6k.lll(kappa, progressive_sieve.r)
    
    g6k.initialize_local(kappa, max(progressive_sieve.r-start_up_n, progressive_sieve.l+1), progressive_sieve.r)
    
    
    progressive_sieve.minl = g6k.l
    progressive_sieve.sat_factor = 1.

    for key in ('g6k', 'verbose','saturation_error'):
       setattr(progressive_sieve, key, locals()[key])

    with g6k.temp_params(reserved_n=progressive_sieve.r-progressive_sieve.l):
        if g6k.params.default_sieve == "gpu":
            progressive_sieve.phase = "up"
        else:
            progressive_sieve.phase = "init"
        wrapped_sieve(progressive_sieve)  # The first initializing Sieve should always be Gauss to avoid rank-loss, sieve process

        progressive_sieve.phase = "up"
        # progressive_sieve Up
        while (g6k.l > progressive_sieve.l):
                
            g6k.extend_left(1)

            if verbose:
                print_pro_randslicer_state(progressive_sieve)
            if not wrapped_sieve(progressive_sieve):
                break
    
    return g6k
            
            