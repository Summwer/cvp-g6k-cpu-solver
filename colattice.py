# from .cvpump_backup import cvpump
# from g6k.utils.stats import dummy_tracer
from copy import deepcopy
# import numpy as np
from fpylll import IntegerMatrix
from randslicer_with_d4f import randslicer_with_d4f
from g6k.siever import Siever
from math import sqrt, log2, log
from strategy_gen.NCSearch import NC_strategy
from fpylll.util import gaussian_heuristic
from time import time

#We should find the close vector from pi_ln(L) to pi_l1(L) where l1<...<ln
def colattice(A, ts, approx_factor, params, sieve_start_index = None, xss2 = None, yl = None, verbose = False, consider_d4f = False, blocksizes =None, optimize_strategy =False, max_sample_times = 1000, min_block = 30, len_bound = 1.):
    """
    Aim to find the closest vector to target vector t.
    
    
    
    Otherwise, sieve_start_index should be d.
    A: lattice matrix
    ts: a set of target vectors 
    approx_factor: approximate factor for approximate (batch)-CVP instance.
    blocksizes: strategy for colattice
    
    
    Case: If we want to find the closest vector to a projected vector, the coefficients of [sieve_start_index, d] is fixed. Then we should set the following parameters: 
    sieve_start_index: should be < d.
    xs: the current coefficients w.r.t lattice basis B.
    yl: the current coefficients w.r.t gso basis B*. 
    len_bound: distance bound for one cvp call, ||w-t||<=len_boung*gh
    """
    
    
    
   
    
    
    g6k = Siever(A,params)
    d = g6k.full_n
    
    
    #Generate blocksize strategy for colattice
    if(blocksizes is None):
        T_sgen = time()
        NC_strategy_class = NC_strategy(min_block = min_block, consider_d4f=consider_d4f, optimize_strategy=optimize_strategy, batch_size = len(ts))
        rr = [g6k.M.get_r(i, i) for i in range(g6k.full_n)]
        gh = sqrt(gaussian_heuristic(rr))
        log_rr = [log(g6k.M.get_r(i, i))/2 for i in range(g6k.full_n)]
        blocksizes, _ = NC_strategy_class.nc_est(log_rr, dist_bound= log(approx_factor*gh))
        T_sgen = time() - T_sgen
        if(max(blocksizes)> 140):
            raise "The sieve dimension exceeds the maximal bound (=140 ), you should do more BKZ preprocess."
    
    
    pts = [[0]*d]*len(ts)
    pws = [[0]*d]*len(ts)
    ws = [[0]*d]*len(ts)
    if sieve_start_index is None:
        sieve_start_index = d
    llb = sieve_start_index
    
    
    xss2 = [[0]*d]*len(ts)
    ees = deepcopy(ts)
    
    if(not consider_d4f):
        f = 0
   
    T_sieve_for_NC = 0
    T_slicer_for_NC = 0
    for i in range(len(blocksizes)-1,-1,-1):
        blocksize = blocksizes[i]
        llb -= blocksize
        if(verbose):
            print("Find close vector on [%d,%d]" %(llb, llb+blocksize))
        
        if(verbose):
            print("Start cvpump-{%d,%d,%d}" %(max(llb,0),blocksize,f))
            
            
        for k in range(len(ts)):
            ts[k] = tuple([ts[k][i] - pts[k][i] + pws[k][i] - ws[k][i] for i in range(d)])
         
        pts, pws, ws, xss, _ , T_sieve , T_slicer, _, gh = randslicer_with_d4f(A,ts, params, kappa=max(llb,0), blocksize = blocksize, f=f, len_bound = len_bound, max_sample_times = max_sample_times)
            
        
        T_sieve_for_NC += T_sieve
        T_slicer_for_NC += T_slicer
            
          
        for k in range(len(ts)):
            ees[k] = tuple([ees[k][i] - ws[k][i] for i in range(d)]) #Final ee = t - sum(w)
            xss2[k] =[xss2[k][i] + xss[k][i] for i in range(d)]
        
        
        
    return blocksizes,ees, xss, T_sgen, T_sieve_for_NC, T_slicer_for_NC, gh

