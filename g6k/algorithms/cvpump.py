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

def print_pro_randslicer_state(cvpump):
    cvpump.minl = min(cvpump.g6k.l, cvpump.minl)
    if cvpump.phase != "down":
        print("\r %3d: ↑%3d      " % (cvpump.g6k.r-cvpump.l, cvpump.g6k.r-cvpump.g6k.l), end=' ')
    else:
        print("\r %3d: ↑%3d ↓%3d " % (cvpump.g6k.r-cvpump.l, cvpump.g6k.r-cvpump.minl, cvpump.g6k.r-cvpump.g6k.l), end=' ')
    sys.stdout.flush()


# Sieve (switching to gauss if needed) and test loop-breaking conditions
# Return true if we should continue
def wrapped_sieve(cvpump):
    if cvpump.phase == "init":
        alg = "gauss"
    else:
        alg = None

    cont = True
    try:
        with cvpump.g6k.temp_params(saturation_ratio=cvpump.g6k.params.saturation_ratio * cvpump.sat_factor):
           
            # Match lifting effort to insertion strategy
            cvpump.g6k(alg=alg, tracer=cvpump.tracer)
            
            
    except SaturationError as e:
        if cvpump.saturation_error == "skip":
            cvpump.down_sieve = False
            logging.info("saturation issue: breaking cvpump.")
            cont = False
        elif cvpump.saturation_error == "weaken":
            logging.info("saturation issue: weakening cvpump.")
            cvpump.sat_factor /= 2.
        elif cvpump.saturation_error == "ignore":
            pass
        else:
            raise e

    # if cvpump.phase == "up" and (cvpump.max_up_time is not None):
    #     if cvpump.max_up_time < time.time() - cvpump.up_time_start:
    #         cont = False

    # if cvpump.goal_r0 is not None:
    #     cvpump.g6k.insert_best_lift(scoring_goal_r0, aux=cvpump)

    # if (cvpump.norm_ee <= cvpump.goal_r0):
    #     cont = False

    return cont


def scoring_goal_r0(i, nlen, olen, aux):
    return i == 0 and nlen < aux.goal_r0

def scoring_down(i, nlen, olen, aux):
    if i < aux.insert_left_bound or nlen >= olen:
        return False
    return log(olen / nlen) - i * log(aux.prefer_left_insert)


def cvpump(g6k, c, tracer, kappa, blocksize,  dim4free,      # Main parameters
         goal_r0=None, start_up_n=30,        
         verbose=False, len_bound = 1                                                          
         ):
    """
    Run the cvpump algorithm.

    :param g6k: The g6k object to work with
    :param c: target vector
    :param tracer: A tracer for g6k
    :param kappa: beginning of the block
    :param blocksize: dimension of the block (r=kappa+blocksize)
    :param dim4free: number of ``dimension for free'' [Ducas, Eurcrypt 2018]: Sieve context [l,r] where l=kappa+dim4free
    :param goal_r0: an extra hook to always insert at position kappa if this goal length can be met
        by a lift.  Quit when this is reached.
    :param verbose: print cvpump steps on the standard output.
    
    return:
       Approximate close vector w to t.
       Differ vector e = t - w.

    """
    cvpump.l = kappa+dim4free  # noqa
    cvpump.r = kappa+blocksize
    cvpump.c = c
    g6k.shrink_db(0)
    g6k.lll(0, cvpump.r)
    # g6k.initialize_local(0, max(g6k.r-start_up_n, cvpump.l+1), l+blocksize)
    g6k.initialize_local(kappa, max(cvpump.r-start_up_n, cvpump.l+1), cvpump.r)
    cvpump.pt = g6k.initialize_target_vector(c)
    # print("pt:",pc)
    
    cvpump.minl = g6k.l
    cvpump.sat_factor = 1.

    for key in ('goal_r0', 'g6k', 'tracer', 'verbose'):
        setattr(cvpump, key, locals()[key])

    with tracer.context(("cvpump", "beta:%d f:%d" % (blocksize, dim4free))):
        with g6k.temp_params(reserved_n=cvpump.r-cvpump.l):
            cvpump.phase = "init"
            wrapped_sieve(cvpump)  # The first initializing Sieve should always be Gauss to avoid rank-loss, sieve process

            cvpump.phase = "up"
            # cvpump Up
            while (g6k.l > cvpump.l):
                if(g6k.n + 1 > g6k.max_sieving_dim):
                    raise RuntimeError("The current sieving context is bigger than maximum supported dimension.")
                g6k.extend_left(1)

                if verbose:
                    print_pro_randslicer_state(cvpump)
                if not wrapped_sieve(cvpump):
                    break
            # print(list(g6k.itervalues()))
            cvpump.pw, cvpump.w, cvpump.x = cvpump.g6k.randslicer(len_bound = len_bound) #pw: projected approximate closest vector; w: approiximate closest vector on full-dim
            # print(cvpump.w)
            # cvpump.ee = [cvpump.c[i] - cvpump.w[i] for i in range(g6k.full_n)]
            # cvpump.norm_ee = norm(cvpump.ee)
            # print(cvpump.ee)

            # if goal_r0 is not None and cvpump.norm_ee <= goal_r0:
            #     return cvpump.w #, cvpump.ee
            
    return cvpump.pt, cvpump.pw, cvpump.w, cvpump.x #, cvpump.ee
