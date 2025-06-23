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
import numpy as np


#the gs-lengths of vector on L[ii:]
def projected_len(ii,vector,M):
    coeffs = list(M.from_canonical(vector))
    coeffs = [0]*ii+coeffs[ii:]
    vec = M.to_canonical(tuple(coeffs))
    return sum([_**2 for _ in vec])

def norm(vec):
    return sum([_**2 for _ in vec])

def print_pro_randslicer_state(bestliftcvp):
    bestliftcvp.minl = min(bestliftcvp.g6k.l, bestliftcvp.minl)
    if bestliftcvp.phase != "down":
        print("\r %3d: ↑%3d      " % (bestliftcvp.g6k.r-bestliftcvp.l, bestliftcvp.g6k.r-bestliftcvp.g6k.l), end=' ')
    else:
        print("\r %3d: ↑%3d ↓%3d " % (bestliftcvp.g6k.r-bestliftcvp.l, bestliftcvp.g6k.r-bestliftcvp.minl, bestliftcvp.g6k.r-bestliftcvp.g6k.l), end=' ')
    sys.stdout.flush()


# Sieve (switching to gauss if needed) and test loop-breaking conditions
# Return true if we should continue
def wrapped_sieve(bestliftcvp):
    if bestliftcvp.phase == "init":
        alg = "gauss"
    else:
        alg = None

    cont = True
    try:
        with bestliftcvp.g6k.temp_params(saturation_ratio=bestliftcvp.g6k.params.saturation_ratio * bestliftcvp.sat_factor):
           
            # Match lifting effort to insertion strategy
            bestliftcvp.g6k(alg=alg, tracer=bestliftcvp.tracer)
            
            
    except SaturationError as e:
        if bestliftcvp.saturation_error == "skip":
            bestliftcvp.down_sieve = False
            logging.info("saturation issue: breaking bestliftcvp.")
            cont = False
        elif bestliftcvp.saturation_error == "weaken":
            logging.info("saturation issue: weakening bestliftcvp.")
            bestliftcvp.sat_factor /= 2.
        elif bestliftcvp.saturation_error == "ignore":
            pass
        else:
            raise e

    # if bestliftcvp.phase == "up" and (bestliftcvp.max_up_time is not None):
    #     if bestliftcvp.max_up_time < time.time() - bestliftcvp.up_time_start:
    #         cont = False

    # if bestliftcvp.goal_r0 is not None:
    #     bestliftcvp.g6k.insert_best_lift(scoring_goal_r0, aux=bestliftcvp)

    # if (bestliftcvp.norm_ee <= bestliftcvp.goal_r0):
    #     cont = False

    return cont


# def scoring_goal_r0(i, nlen, olen, aux):
#     return i == 0 and nlen < aux.goal_r0

# def scoring_down(i, nlen, olen, aux):
#     if i < aux.insert_left_bound or nlen >= olen:
#         return False
#     return log(olen / nlen) - i * log(aux.prefer_left_insert)



def to_canonical(g6k,yr,kappa):
    d = len(yr)
    # dim = d - start
    # 1. triangular system solving
    for i in range(d-1,kappa-1,-1):
        for j in range(i+1,d):    
            yr[i] -= g6k.M.get_mu(j, i) * yr[j]
            # if(i==5):
            #     print(g6k.M.get_mu(j, i),end=",")
            # print("mu:", g6k.M.get_mu(start + j, start + i),g6k.M.get_mu(start + i, start + j))
    # print("v_yr:",yr)
    # 2. multiply by B
    return g6k.M.B.multiply_left(tuple(yr))



def bestliftcvp(g6k, g6k2, insert_left_bound, tracer, kappa, blocksize,  dim4free,      # Main parameters
         goal_r0=None, start_up_n=30,        
         verbose=False, len_bound = 1                                                          
         ):
    """
    Run the bestliftcvp algorithm.

    :param g6k: The g6k object to work with for slicer
    :param g6k2: The g6k2 object to work with for storing sieve information
    :param L:  best lift List, contains multiple target vectors under same lattice basis.
    :param tracer: A tracer for g6k
    :param kappa: beginning of the block
    :param blocksize: dimension of the block (r=kappa+blocksize)
    :param dim4free: number of ``dimension for free'' [Ducas, Eurcrypt 2018]: Sieve context [l,r] where l=kappa+dim4free
    :param goal_r0: an extra hook to always insert at position kappa if this goal length can be met
        by a lift.  Quit when this is reached.
    :param verbose: print bestliftcvp steps on the standard output.
    
    return:
       Approximate close vector w to t.
       Differ vector e = t - w.

    """
    bestliftcvp.l = kappa+dim4free  # noqa
    bestliftcvp.r = kappa+blocksize
    
    g6k.shrink_db(0)
    # g6k.lll(0, bestliftcvp.r)
    # g6k.initialize_local(0, max(g6k.r-start_up_n, bestliftcvp.l+1), l+blocksize)
    g6k.initialize_local(kappa, max(bestliftcvp.r-start_up_n, bestliftcvp.l+1), bestliftcvp.r)
    

    # print("pt:",pc)
    
    bestliftcvp.minl = g6k.l
    bestliftcvp.sat_factor = 1.

    for key in ('goal_r0', 'g6k', 'tracer', 'verbose'):
        setattr(bestliftcvp, key, locals()[key])

    with tracer.context(("bestliftcvp", "beta:%d f:%d" % (blocksize, dim4free))):
        with g6k.temp_params(reserved_n=bestliftcvp.r-bestliftcvp.l):
            bestliftcvp.phase = "init"
            wrapped_sieve(bestliftcvp)  # The first initializing Sieve should always be Gauss to avoid rank-loss, sieve process

            bestliftcvp.phase = "up"
            # bestliftcvp Up
            while (g6k.l > bestliftcvp.l):
                if(g6k.n + 1 > g6k.max_sieving_dim):
                    raise RuntimeError("The current sieving context is bigger than maximum supported dimension.")
                g6k.extend_left(1)

                if verbose:
                    print_pro_randslicer_state(bestliftcvp)
                if not wrapped_sieve(bestliftcvp):
                    break
            # print(list(g6k.itervalues()))
            
            
            
            minnorm = None
            ii = None
            best_v = None
            minee = None
            # j = 0
            
            L = [ (index, nlen, v, yr) for (index, nlen, v, yr) in g6k2.best_lifts() if index >=insert_left_bound]
            Ll = g6k2.l
            Lr = g6k2.r
            
            #db = list(g6k.itervalues())
            for (index, nlen, v , yr) in L:

                w  = g6k.M.B.multiply_left(tuple([0]*Ll + list(v)[Ll:Lr]))
                yr = [0]*Ll + list(yr)[Ll:Lr]
                pw = g6k.M.to_canonical(yr)          
                t= tuple([pw[i] - w[i] for i in range(Lr)])

                # print([round(_,2) for _ in g6k.M.from_canonical(tuple(t))])
                
                g6k.initialize_target_vector(t,0)
                
                #print(t)
                # print(w)

                bestliftcvp.pw, bestliftcvp.w, bestliftcvp.x = bestliftcvp.g6k.randslicer(len_bound = len_bound) 
                
                ee = tuple([w[i] + bestliftcvp.w[i] for i in range(Lr)])
                
                # print("inserted vector:", ee)
                # print("inserted vector coeff:",[round(_,2) for _ in list(np.linalg.solve((np.array(list(g6k.M.B))).T,np.array(ee)))] )
            
                norm_ee = sum([_**2 for _ in ee]) 
                
                
                #if(index == 0):
                
                
                if((minnorm is None or minnorm > norm_ee) and norm_ee > 0):
                    # print(index,nlen,norm_ee,goal_r0)
                    minnorm = norm_ee
                    minee = ee
                    best_v = [0]*Lr
                    for i in range(Lr):
                        if(i<Ll):
                            best_v[i] = bestliftcvp.x[i]
                        else:
                            best_v[i] = (bestliftcvp.x[i]+v[i])
                        best_v = np.array(best_v)
                    # print(best_v)
                    ii = index
                    for i in range(index):
                        if(projected_len(ii,ee,g6k.M)<g6k.M.get_r(i,i)):
                            ii = i
                            break
                    if(goal_r0 is not None and norm_ee<=goal_r0):
                        break
                # break
    # print("ii:", ii)
    # print("best_v:", best_v)            
    print("minnorm:", minnorm)
    # print("ee:", tuple(np.dot(best_v,Ba)))
    return ii, best_v,minnorm,minee #, bestliftcvp.ee
