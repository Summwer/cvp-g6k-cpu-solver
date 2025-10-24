import time
from g6k.algorithms.progressive_sieve_preprocess import progressive_sieve
# from g6k.algorithms.cvpump import cvpump
from math import sqrt
from fpylll import *
from fpylll.util import gaussian_heuristic
from fpylll import *
from g6k.siever import Siever
from g6k.utils.stats import dummy_tracer
from multiprocessing import Pool
from functools import partial
from g6k.siever_params import SieverParams
from copy import deepcopy
from time import sleep




def dot_product(a,b):
    return sum([a[i]*b[i] for i in range(len(a))])

def norm(a):
    return sqrt(dot_product(a,a))

def cvp_enum(A, kappa, blocksize, t):

    n = len(t)
    
    
    
    M = GSO.Mat(A, float_type = 'ld')
    M.update_gso()
    yl = M.from_canonical(t)
   
    pt_yl = [0]*n
    for i in range(kappa, kappa+blocksize):
        pt_yl[i] = yl[i]
    pt = M.to_canonical(tuple(pt_yl))
    
    rr = [M.get_r(i, i) for i in range(n)]
    
    
    enum_obj = Enumeration(M, nr_solutions=1)
    res = enum_obj.enumerate(kappa, kappa+blocksize, 2*gaussian_heuristic(rr[kappa:kappa+blocksize]), 0, target=pt_yl)
    v = res[0][1]

    mat_x = IntegerMatrix(1,n)
    for i in range(kappa, kappa+blocksize):
        mat_x[0,i] = int(v[i-kappa])


    w =  mat_x * M.B  
    x = list(mat_x[0])
    w =  list(w[0])
    pw_yl = [0]*(kappa) + list(M.from_canonical(tuple(w)))[kappa:kappa+blocksize] + [0]*(n-blocksize-kappa)
    pw = M.to_canonical(pw_yl)
  
    return pt,pw, w, x


def randslicer_with_d4f(A,ts, params, kappa=0, blocksize = None, f=0, len_bound = 1., max_sample_times = 1000, verbose = False):
    """
    Use randomized slicer to solve (approximate) CVP and bath-CVP problem.
    A: lattice basis.
    ts: list of target vectors.
    kappa: left bound
    blocksize: cvp range
    f: dimension-for-free value. 
    """
    # params = SieverParams(threads=threads, 
                                # default_sieve=default_sieve, reserved_db_size= max((2 * sqrt(4/3.) ** len(ts[0])), 2000))
    #
    
    
    if(blocksize is None):
        blocksize = g6k.full_n
    params.otf_lift = False
    g6k = Siever(A,params)
    rr = [g6k.M.get_r(i, i) for i in range(g6k.full_n)]
    gh = sqrt(gaussian_heuristic(rr[kappa:kappa+blocksize]))
    
    
    
    
    #w = tuple([0]*A.ncols) 
    sample_times = [1]*len(ts)
    T_sieve = 0
    
    
    pts = []
    ws = []
    pws = []
    xs = []
    
    if(blocksize <= 45):
        #CVP enumeration 
        T0 =time.time()
        # print(kappa, blocksize)
        fixed_params = (A,kappa, blocksize)
        processor = partial(cvp_enum, *fixed_params)
        with Pool(processes=g6k.params.threads) as pool:
            for pt,pw, w, x in pool.imap(processor, ts):
                pts.append(pt)
                ws.append(w)
                pws.append(pw)
                xs.append(x)
        # pool.close()
        T_slicer = time.time() - T0
        
        
    else:
    # if (True):
        T_sieve = time.time()
        progressive_sieve(g6k, kappa, blocksize,  f, verbose = verbose)
        #cvpump(g6k,dummy_tracer,max(kappa,0),blocksize, f,verbose=True)
        T_sieve = time.time() - T_sieve
       
        
        T_bucket = time.time()
        g6k.cdb_bucket_process()
        T_bucket = time.time() - T_bucket
        # print(g6k.db_size())
        
        
        
            
        
       
        T_slicer = time.time()
        pts = g6k.initialize_target_vectors(ts)
        (pws,ws,xs),sample_times  = g6k.randslicer(max_sample_times = max_sample_times,len_bound = len_bound)    
        T_slicer = time.time()-T_slicer
       
                
        # print("pws = ", pws)
        # print("ws = ", ws)
        # print("xs = ", xs)
        # print("sample times = ", sample_times)
        
        # for t in ts:
        #     pt = g6k.initialize_target_vectors(ts)
        #     pw,w,x,sample_times  = g6k.randslicer(max_sample_times = max_sample_times,len_bound = len_bound)
        #     # pw, w, x, sample_times = cvpump.g6k.randslicer(len_bound = len_bound)
        #     # pw, w, x = g6k.get_cv()
        #     # print("sample_times:", sample_times)
        #     pts.append(pt) 
        #     pws.append(pw)
        #     ws.append(w)
        #     xs.append(x)
            # print("pw = ", pw)
            # print("pt = ", pt)
        
        
         #need to modify
    return pts, pws, ws,xs, sample_times, T_sieve,  T_slicer, g6k.db_size(), gh
