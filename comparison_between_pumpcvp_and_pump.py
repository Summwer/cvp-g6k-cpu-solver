from g6k.utils.util import load_svpchallenge_and_randomize
from g6k.algorithms.pumpcvp import pumpcvp
from g6k.algorithms.pump import pump
from g6k.siever import Siever
from g6k.utils.stats import SieveTreeTracer, dummy_tracer
from fpylll.tools.quality import basis_quality
from fpylll.util import gaussian_heuristic
import time

def svpchal_pumpcvp_test(n, dsvp,verbose=True,tracer=dummy_tracer):
    A = None #load_svp_midmat(n)
    if(A is None):
        A, _ = load_svpchallenge_and_randomize(n, s=0, seed=0)
    if verbose:
        print(("Loaded challenge dim %d" % n))
        
    g6k = Siever(A,None)
    g6k.lll(0,g6k.full_n)
    gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(n)])
    goal_r0 = (1.05 ** 2) * gh
    print("gh = %d, goal_r0 = %.3f" %(gh,goal_r0))
    slope = basis_quality(g6k.M)["/"]
    print("Intial Slope = %.5f\n" % slope)
    f = n - dsvp
    print("pumpcvp-{0,%d,%d}, sieve_dim = %d" %(n,f,dsvp))
    print("rr[0]=", g6k.M.get_r(0,0))
    T0 = time.time()
    pumpcvp(g6k,tracer, 0, n, f, verbose=True, down_sieve = True, min_cvp_dim = 50, goal_r0 = goal_r0)
    g6k.lll(0,g6k.full_n)
    g6k.M.update_gso()
    norms = []
    for i in range(g6k.full_n):
        norm_i  = 0
        for j in range(g6k.full_n):
            norm_i += g6k.M.B[i][j]*g6k.M.B[i][j]
        norms.append(norm_i)
    print(min(norms))
    print("rr[0]=", g6k.M.get_r(0,0))
    print("walltime cost for pumpcvp-{0,%d,%d}: %.2f sec\n" %(n,f,time.time()-T0))


def svpchal_pump_test(n, dsvp,verbose=True,tracer=dummy_tracer):
    A = None #load_svp_midmat(n)
    if(A is None):
        A, _ = load_svpchallenge_and_randomize(n, s=0, seed=0)
    if verbose:
        print(("Loaded challenge dim %d" % n))
        
    f = n - dsvp
    g6k = Siever(A,None)
    gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(n)])
    goal_r0 = (1.05 ** 2) * gh
    print("gh = %d, goal_r0 = %.3f" %(gh,goal_r0))
    slope = basis_quality(g6k.M)["/"]
    print("Intial Slope = %.5f\n" % slope)
    g6k.lll(0,g6k.full_n)
    print("Pump-{0,%d,%d}, sieve_dim = %d" %(n,f,dsvp))
    T0 = time.time()
    pump(g6k,tracer, 0, n, f, verbose=True,  down_sieve = True, goal_r0=goal_r0)
    g6k.lll(0,g6k.full_n)
    print("rr[0]=", g6k.M.get_r(0,0))
    print("walltime cost for pump-{0,%d,%d}: %.2f sec\n" %(n,f,time.time()-T0))
    
    
    
    
n = 140
dsvp = n//2
svpchal_pumpcvp_test(n, dsvp-10)
svpchal_pump_test(n, dsvp)