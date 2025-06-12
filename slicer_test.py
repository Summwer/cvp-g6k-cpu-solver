from fpylll.util import gaussian_heuristic
from fpylll import *
from g6k.siever import Siever
from g6k.utils.stats import dummy_tracer
from g6k.siever_params import SieverParams
from g6k.algorithms.pro_randslicer import pro_randslicer
from math import sqrt
from g6k.utils.util import load_cvp_instance
from math import ceil, exp
import time
from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
from fpylll.tools.quality import basis_quality

def dot_product(a,b):
    return sum([a[i]*b[i] for i in range(len(a))])

def norm(a):
    return sqrt(dot_product(a,a))


#Use the simulator like BKZ 2.0
def draw_cvp_bound_simulation():
    return


def cvp_test(A,t, params, len_bound = 1., max_sample_times = 1000):
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
        # params = SieverParams(threads = 1,saturation_ratio = 1.)
        
        g6k = Siever(A,params)
        f = 0
        # print(g6k.M.B.nrows,g6k.M.B.ncols)
        close_vector,_,sample_times, T_sieve, T_slicer, db_size = pro_randslicer(g6k,t,dummy_tracer,f,verbose=False, len_bound = len_bound, max_sample_times = max_sample_times)
    return close_vector, sample_times, T_sieve, T_slicer, db_size





#set alpha = sqrt(4/3), the success probability is: exp(-0.149780140929454 * beta)


print("{0: <10} | {1: <10} | {2: <15} | {3: <30} | {4: <15} | {5: <15} | {6: <15} | {7: <15} | {8: <15} | {9: <15}".format("dim", "index", "sample times", "estimated sample times", "T_pump (sec)", "T_slicer (sec)", "dt", "gh", "db_size", "satisfied vectors"))
tours = 10



for n in range(60, 100, 5):
    for index in range(1,tours+1):
        A, t,len_bound = load_cvp_instance(n, approx_factor = 0.99)
        
        
        sp = SieverParams()
        
        params = SieverParams(threads = 1 , default_sieve = "bdgl2", reserved_db_size = sp["db_size_factor"] * sp["db_size_base"] ** n)
        
        g6k = Siever(A,params)
        
            # slope = basis_quality(g6k.M)["/"]
            # fmt = "slope: %.5f, walltime: %.3f sec"
            # print(fmt % (slope, time.time() - T0))
        rr = [g6k.M.get_r(i, i) for i in range(n)]

        
        max_sample_times = 140
        w, sample_times,T_pump, T_slicer, db_size = cvp_test(g6k.M.B,t, params, len_bound = len_bound, max_sample_times = max_sample_times)
        # max_sample_times = ceil((16/13.)**(n//2.))
         #ceil(exp(0.149780140929454 * n))
        
        gh = sqrt(gaussian_heuristic(rr))
        dt = sqrt(sum([(w[i] - t[i])**2 for i in range(len(t))]))
        # simDist = DistEstColattice([log(_)/2. for _ in rr[:n]], [n])
    
        print("{0: <10} | {1:<10} | {2: <15} | {3: <30} | {4: <15} | {5: <15} | {6: <15} | {7: <15} | {8: <15} | {9: <15}".format(n,index, sample_times, max_sample_times, round(T_pump,4), round(T_slicer,4), round(dt,3), round(gh,3), db_size, int(.5 * params.saturation_ratio * params.db_size_base ** n )))

