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
from g6k.siever import SaturationError, Siever
import logging
import numpy as np
from g6k.algorithms.bestliftcvp import bestliftcvp
from g6k.utils.stats import dummy_tracer
from copy import deepcopy

def print_pumpcvp_state(pumpcvp):
    pumpcvp.minl = min(pumpcvp.g6k.l, pumpcvp.minl)
    if pumpcvp.phase != "down":
        print("\r %3d: ↑%3d      " % (pumpcvp.r-pumpcvp.l, pumpcvp.g6k.r-pumpcvp.g6k.l), end=' ')
    else:
        print("\r %3d: ↑%3d ↓%3d " % (pumpcvp.r-pumpcvp.l, pumpcvp.r-pumpcvp.minl, pumpcvp.r-pumpcvp.g6k.l), end=' ')
    sys.stdout.flush()


# Sieve (switching to gauss if needed) and test loop-breaking conditions
# Return true if we should continue
def wrapped_sieve(pumpcvp):
    if pumpcvp.phase == "init":
        alg = "gauss"
    else:
        alg = None

    cont = True
    try:
        with pumpcvp.g6k.temp_params(saturation_ratio=pumpcvp.g6k.params.saturation_ratio * pumpcvp.sat_factor):
            # Match lifting effort to insertion strategy
            pumpcvp.g6k(alg=alg, tracer=pumpcvp.tracer)
            
    except SaturationError as e:
        if pumpcvp.saturation_error == "skip":
            pumpcvp.down_sieve = False
            logging.info("saturation issue: breaking pumpcvp.")
            cont = False
        elif pumpcvp.saturation_error == "weaken":
            logging.info("saturation issue: weakening pumpcvp.")
            pumpcvp.sat_factor /= 2.
        elif pumpcvp.saturation_error == "ignore":
            pass
        else:
            raise e

    if pumpcvp.phase == "up" and (pumpcvp.max_up_time is not None):
        if pumpcvp.max_up_time < time.time() - pumpcvp.up_time_start:
            cont = False

    if pumpcvp.goal_r0 is not None:
        pumpcvp.g6k.insert_best_lift(scoring_goal_r0, aux=pumpcvp)

        if (pumpcvp.g6k.M.get_r(pumpcvp.kappa, pumpcvp.kappa) <= pumpcvp.goal_r0):
            cont = False

    return cont


def scoring_goal_r0(i, nlen, olen, aux):
    return i == aux.kappa and nlen < aux.goal_r0


def scoring_down(i, nlen, olen, aux):
    if i < aux.insert_left_bound or nlen >= olen:
        return False
    return log(olen / nlen) - i * log(aux.prefer_left_insert)


def pumpcvp(g6k, tracer, kappa, blocksize, dim4free, down_sieve=False,                                 # Main parameters
         goal_r0=None, max_up_time=None, down_stop=None, start_up_n=30, saturation_error="weaken",  # Flow control of the pumpcvp
         increasing_insert_index=True, prefer_left_insert=1.04,                                     # Insertion policy
         verbose=False, min_cvp_dim = 70                                                                         # Misc
         ):
    """
    Run the pumpcvp algorithm.

    :param g6k: The g6k object to work with
    :param tracer: A tracer for g6k
    :param kappa: beginning of the block
    :param blocksize: dimension of the block (r=kappa+blocksize)
    :param dim4free: number of ``dimension for free'' [Ducas, Eurcrypt 2018]: Sieve context [l,r] where l=kappa+dim4free
    :param down_sieve: re-sieve after each insert during the pumpcvp-down phase.  (stronger reduction,
        slower running time)
    :param goal_r0: an extra hook to always insert at position kappa if this goal length can be met
        by a lift.  Quit when this is reached.
    :param max_up_time: For balancing BKZ time with SVP call time in LWE.  Stop pumpcvping up when this
        time has elapsed.
    :param down_stop: stop inserts during pumpcvping down after index kappa+down_stop (to control
        overheads of insert in large dimensional lattices)
    :param start_up_n: Initial sieve-context dimension for pumpcvping up (starting at 1 incurs useless overheads)
    :param saturation_error: determines the behavior of pumpcvp when encountering saturation issue {"weaken",
        "skip", "ignore", "forward"}
    :param increasing_insert_index: During pumpcvp-down, always insert on the right side of the previous insertion.
    :param prefer_left_insert: Parameter theta from the paper (Sec 4.4) for scoring insertion candidates.
    :param verbose: print pumpcvp steps on the standard output.

    """
    pumpcvp.r = kappa+blocksize
    pumpcvp.l = kappa+dim4free  # noqa

    max_cvp_dim = max(kappa+blocksize,blocksize - dim4free)

    g6k.shrink_db(0)
    g6k.lll(kappa, pumpcvp.r)
    g6k.initialize_local(kappa, max(pumpcvp.r-start_up_n, pumpcvp.l+1), pumpcvp.r)

    pumpcvp.sat_factor = 1.
    pumpcvp.up_time_start = time.time()
    pumpcvp.insert_left_bound = kappa
    pumpcvp.minl = g6k.l

    for key in ('kappa', 'down_sieve', 'goal_r0', 'g6k', 'tracer',
                'max_up_time', 'saturation_error', 'verbose', 'prefer_left_insert'):
        setattr(pumpcvp, key, locals()[key])

    if down_stop is None:
        down_stop = dim4free

    with tracer.context(("pumpcvp", "beta:%d f:%d" % (blocksize, dim4free))):
        with g6k.temp_params(reserved_n=pumpcvp.r-pumpcvp.l):
            pumpcvp.phase = "init"
            wrapped_sieve(pumpcvp)  # The first initializing Sieve should always be Gauss to avoid rank-loss

            pumpcvp.phase = "up"
            # pumpcvp Up
            while (g6k.l > pumpcvp.l):
                if(g6k.n + 1 > g6k.max_sieving_dim):
                    raise RuntimeError("The current sieving context is bigger than maximum supported dimension.")
                g6k.extend_left(1)

                if verbose:
                    print_pumpcvp_state(pumpcvp)
                if not wrapped_sieve(pumpcvp):
                    break

            if goal_r0 is not None and (g6k.M.get_r(kappa, kappa) <= goal_r0):
                return

            # pumpcvp Down
            pumpcvp.phase = "down"



            tmpg6k = Siever(g6k.M.B,None)
            while (g6k.n > 3 and pumpcvp.insert_left_bound <= kappa+down_stop):
                # (try to) Insert      
                ii = None
                # if(g6k.l < max_cvp_dim):
                if(g6k.l == pumpcvp.l):
                # if(True):
                    L = [ (index, nlen, v) for (index, nlen, v) in g6k.best_lifts() if index >=pumpcvp.insert_left_bound]
                    print("best lift list size: ", len(L))
                    if(len(L)!=0):
                        ii,best_v,minnorm,minee = bestliftcvp(tmpg6k,L,g6k.l,g6k.r,dummy_tracer,0, max(g6k.l,min_cvp_dim), 0,goal_r0=goal_r0,verbose=verbose)#, len_bound = 0.5)            

                        print(ii,pumpcvp.insert_left_bound)
                        
                        if(ii is not None):
                            g6k.insert_cvp(ii, best_v)
                            # g6k.insert(ii, best_v)
                            # break
                        if(minnorm<goal_r0):
                            print("solution:",minee)
                            return 
                else:
                    ii = g6k.insert_best_lift(scoring_down, aux=pumpcvp)
                
                print("rr[0]:",g6k.M.get_r(0,0))

                if ii is not None and increasing_insert_index:
                    pumpcvp.insert_left_bound = ii + 1
                else:
                    g6k.shrink_left(1)
                
                # print(g6k.M.get_r(kappa, kappa))
                if goal_r0 is not None and (g6k.M.get_r(kappa, kappa) <= goal_r0):
                    break

                # Sieve (or Shrink db)
                if verbose:
                    print_pumpcvp_state(pumpcvp)
                if not pumpcvp.down_sieve:
                    g6k.resize_db(max(500, g6k.db_size() / g6k.params.db_size_base))
                elif not wrapped_sieve(pumpcvp):
                    break
