from fpylll.util import gaussian_heuristic
from fpylll import *
from g6k.siever import Siever
from g6k.utils.stats import dummy_tracer
from g6k.siever_params import SieverParams
from g6k.algorithms.pro_randslicer import pro_randslicer
from math import sqrt
from g6k.utils.util import load_cvp_instance
from copy import deepcopy
from DistEstColattice import DistEstDistEstColattice
from math import log, ceil
from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
import time

def dot_product(a,b):
    return sum([a[i]*b[i] for i in range(len(a))])

def norm(a):
    return sqrt(dot_product(a,a))


#Use the simulator like BKZ 2.0
def draw_cvp_bound_simulation():
    return


def cvp_test(A,t):
    close_vector = tuple([0]*A.ncols) 
    sample_times = 1
    T_sieve = 0
    db_size = 0
    if(A.nrows == 1):
        T0 =time.time()
        #babai
        c = round(dot_product(t,A[0])/ dot_product(A[0],A[0]))
        close_vector = tuple([c*A[0][i] for i in range(A.ncols)]) 
        T_slicer = time.time() - T0
    elif(A.nrows <= 30):
        #CVP enumeration 
        A = IntegerMatrix.from_matrix(A, int_type="mpz")
        T0 =time.time()
        close_vector = CVP.closest_vector(A,t)
        T_slicer = time.time() - T0
    else:
        #randomlized slicer
        params = SieverParams(threads = 1)
        g6k = Siever(A,params)
        f = 0
        # print(g6k.M.B.nrows,g6k.M.B.ncols)
        close_vector,_,sample_times, T_sieve, T_slicer, db_size = pro_randslicer(g6k,t,dummy_tracer,f,verbose=False)
    return close_vector, sample_times, T_sieve, T_slicer, db_size




print("{0: <10} | {1: <10} | {2: <15} | {3: <30} | {4: <15} | {5: <15} | {6: <15} | {7: <15} | {8: <15}".format("dim", "index", "sample times", "estimated sample times", "T_pump (sec)", "T_slicer (sec)", "dt", "gh", "db_size"))
tours = 10
for n in range(50, 100, 2):
    for index in range(1,tours+1):
        A, t = load_cvp_instance(n)
        A = LLL.reduction(A)
        
        
        g6k = Siever(A,None)
        for blocksize in range(10, 30):
            bkz = BKZReduction(g6k.M)
            par = fplll_bkz.Param(blocksize,
                                        strategies=fplll_bkz.DEFAULT_STRATEGY,
                                        max_loops=1)
            bkz(par)
        
        
        rr = [g6k.M.get_r(i, i) for i in range(n)]


        t_yl = g6k.M.from_canonical(t)
        pt = g6k.M.to_canonical(tuple(list(t_yl)))
        
        w, sample_times,T_pump, T_slicer, db_size = cvp_test(A,t)
        max_sample_times = ceil((16/13.)**(n//2.))
        
        
        gh = sqrt(gaussian_heuristic(rr))
        dt = sqrt(sum([(w[i] - pt[i])**2 for i in range(len(t))]))
        # simDist = DistEstDistEstColattice([log(_)/2. for _ in rr[:n]], [n])
    
        print("{0: <10} | {1:<10} | {2: <15} | {3: <30} | {4: <15} | {5: <15} | {6: <15} | {7: <15} | {8: <15}".format(n,index, sample_times, max_sample_times, round(T_pump,4), round(T_slicer,4), round(dt,3), round(gh,3), db_size))

