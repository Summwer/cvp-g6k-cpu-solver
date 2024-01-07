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


# def initialize_target_vector(long* target_vector, fplll::MatGSO<SZT, SFT> M, vector<SFT> &yl){
# //     vector<SFT> z = vector<SFT>(full_n,0.); //use yl to store the <t,bi*>/<bi*,bi*>
# //     //convert the type from float to SFT
# //     for(unsigned int i = 0; i < full_n; i++){
# //       z[i] = target_vector[i]; 
# //     }
# //     M.from_canonical(yl,z);  //y is the coeefficient of target_vector wrt B*.
# // }

def norm(vec):
    return sum([_**2 for _ in vec])

def print_pro_randslicer_state(pro_randslicer):
    pro_randslicer.minl = min(pro_randslicer.g6k.l, pro_randslicer.minl)
    if pro_randslicer.phase != "down":
        print("\r %3d: ↑%3d      " % (pro_randslicer.g6k.full_n-pro_randslicer.l, pro_randslicer.g6k.r-pro_randslicer.g6k.l), end=' ')
    else:
        print("\r %3d: ↑%3d ↓%3d " % (pro_randslicer.g6k.full_n-pro_randslicer.l, pro_randslicer.g6k.full_n-pro_randslicer.minl, pro_randslicer.g6k.full_n-pro_randslicer.g6k.l), end=' ')
    sys.stdout.flush()


# Sieve (switching to gauss if needed) and test loop-breaking conditions
# Return true if we should continue
def wrapped_sieve(pro_randslicer):
    if pro_randslicer.phase == "init":
        alg = "gauss"
    else:
        alg = None

    cont = True
    try:
        with pro_randslicer.g6k.temp_params(saturation_ratio=pro_randslicer.g6k.params.saturation_ratio * pro_randslicer.sat_factor):
           
            # Match lifting effort to insertion strategy
            pro_randslicer.g6k(alg=alg, tracer=pro_randslicer.tracer)
            
            
    except SaturationError as e:
        if pro_randslicer.saturation_error == "skip":
            pro_randslicer.down_sieve = False
            logging.info("saturation issue: breaking pro_randslicer.")
            cont = False
        elif pro_randslicer.saturation_error == "weaken":
            logging.info("saturation issue: weakening pro_randslicer.")
            pro_randslicer.sat_factor /= 2.
        elif pro_randslicer.saturation_error == "ignore":
            pass
        else:
            raise e

    # if pro_randslicer.phase == "up" and (pro_randslicer.max_up_time is not None):
    #     if pro_randslicer.max_up_time < time.time() - pro_randslicer.up_time_start:
    #         cont = False

    # if pro_randslicer.goal_r0 is not None:
    #     pro_randslicer.g6k.insert_best_lift(scoring_goal_r0, aux=pro_randslicer)

    # if (pro_randslicer.norm_ee <= pro_randslicer.goal_r0):
    #     cont = False

    return cont


def scoring_goal_r0(i, nlen, olen, aux):
    return i == 0 and nlen < aux.goal_r0

def scoring_down(i, nlen, olen, aux):
    if i < aux.insert_left_bound or nlen >= olen:
        return False
    return log(olen / nlen) - i * log(aux.prefer_left_insert)


def pro_randslicer(g6k, c, tracer, dim4free,        # Main parameters
         goal_r0=None, start_up_n=30,        
         verbose=False,                                                               
         ):
    """
    Run the pro_randslicer algorithm.

    :param g6k: The g6k object to work with
    :param c: target vector
    :param tracer: A tracer for g6k
    :param dim4free: number of ``dimension for free'' [Ducas, Eurcrypt 2018]: Sieve context [l,r] where l=dim4free
    :param goal_r0: an extra hook to always insert at position kappa if this goal length can be met
        by a lift.  Quit when this is reached.
    :param verbose: print pro_randslicer steps on the standard output.
    
    return:
       Approximate close vector w to t.
       Differ vector e = t - w.

    """
    pro_randslicer.l = dim4free  # noqa
    pro_randslicer.c = c
    g6k.shrink_db(0)
    g6k.lll(0, g6k.full_n)
    g6k.initialize_local(0, max(g6k.full_n-start_up_n, pro_randslicer.l+1), g6k.full_n)
    g6k.initialize_target_vector(c)
    
    
    
    pro_randslicer.minl = g6k.l
    pro_randslicer.sat_factor = 1.

    for key in ('goal_r0', 'g6k', 'tracer', 'verbose'):
        setattr(pro_randslicer, key, locals()[key])

    with tracer.context(("pro_randslicer", "d:%d f:%d" % (g6k.full_n, dim4free))):
        with g6k.temp_params(reserved_n=g6k.full_n-pro_randslicer.l):
            pro_randslicer.phase = "init"
            wrapped_sieve(pro_randslicer)  # The first initializing Sieve should always be Gauss to avoid rank-loss, sieve process

            pro_randslicer.phase = "up"
            # pro_randslicer Up
            while (g6k.l > pro_randslicer.l):
                if(g6k.n + 1 > g6k.max_sieving_dim):
                    raise RuntimeError("The current sieving context is bigger than maximum supported dimension.")
                g6k.extend_left(1)

                if verbose:
                    print_pro_randslicer_state(pro_randslicer)
                if not wrapped_sieve(pro_randslicer):
                    break
            # print(list(g6k.itervalues()))
            pro_randslicer.w = pro_randslicer.g6k.randslicer()
            pro_randslicer.ee = [pro_randslicer.c[i] - pro_randslicer.w[i] for i in range(pro_randslicer.g6k.full_n)]
            pro_randslicer.norm_ee = norm(pro_randslicer.ee)
            # print(pro_randslicer.ee)

            if goal_r0 is not None and pro_randslicer.norm_ee <= goal_r0:
                return pro_randslicer.w, pro_randslicer.ee
            
    return pro_randslicer.w, pro_randslicer.ee
