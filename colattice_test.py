
from g6k.siever_params import SieverParams
from math import sqrt
from g6k.utils.util import load_cvp_instance
# from full_cvp_solver import full_cvp_solver
import argparse
from colattice import colattice


def colattice_test(min_dim, max_dim, step, tours, approx_factor, 
               threads, default_sieve, max_sample_times, batch_size, len_bound, optimize_strategy ):
    """Run Colattice solver tests with given parameters"""
    print("{0: <10} | {1: <50} | {2: <15} | {3: <15} | {4: <20} | {5: <15} | {6: <15} | {7: <15}".format(
        "dim", "strategy", "index", "T_sgen (sec)",
        "T_sieve (sec)", "T_colattice (sec)", "len_bound", "required len_bound"))
    
    for n in range(min_dim, max(max_dim+1, min_dim+1), step):
        for index in range(1, tours+1):
            A, ts  = load_cvp_instance(n, batch_size = batch_size)
            
            #sp = SieverParams()
            
            params = SieverParams(threads=threads, 
                                default_sieve=default_sieve)#, saturation_ratio = 2.)
            
            max_sample_times = max_sample_times
            strategy, ees, _, T_sgen, T_pump, T_colattice, gh = colattice(A, ts, approx_factor, params, len_bound = len_bound, optimize_strategy = optimize_strategy )
            
            dt_min = None
            dt_max = None
            for ee in ees:
                dt = sqrt(sum([ee[i]**2 for i in range(len(ee))]))
                if(dt_min is None or dt < dt_min):
                    dt_min = dt 
                if(dt_max is None or dt > dt_max):
                    dt_max = dt 
            
            if(len(ts) == 1):
                actual_len_bound = round(dt_min/gh,3)
            else:
                actual_len_bound = str( round(dt_min/gh,3))+"~" + str( round(dt_max/gh,3))
            print("{0: <10} | {1:<50} | {2: <15} | {3: <15} | {4: <20} | {5: <15} | {6: <15} | {7: <15}".format(
                    n, str(strategy), index, round(T_sgen, 4), 
                    round(T_pump,4), round(T_colattice,4), 
                    actual_len_bound, round(approx_factor,3)))
            

#Use the simulator like BKZ 2.0
def draw_cvp_bound_simulation():
    return






if __name__ == '__main__':
    description = "Run Colattice solver tests with configurable parameters"
    
    # Setup argument parser
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--min-dim', type=int, default=60,
                       help='Minimum dimension (default: 60)')
    parser.add_argument('--max-dim', type=int, default=100,
                       help='Maximum dimension (default: 100)')
    parser.add_argument('--step', type=int, default=5,
                       help='Dimension step size (default: 5)')
    parser.add_argument('--tours', type=int, default=10,
                       help='Number of tours per dimension (default: 10)')
    parser.add_argument('--approx-factor', type=float, default=0.99,
                       help='Approximation factor (default: 0.99)')
    parser.add_argument('--threads', type=int, default=1,
                       help='Number of threads (default: 1)')
    parser.add_argument('--default-sieve', type=str, default="bdgl2",
                       help='Default sieve algorithm (default: bdgl2)')
    parser.add_argument('--max-sample-times', type=int, default=1000,
                       help='Maximum sample times (default: 1000)')
    parser.add_argument('--batch-size', type=int, default=1,
                       help='number of target vector (default: 1)')
    parser.add_argument('--len-bound', type=float, default=1.,
                       help='the square distance bound w.r.t gh for one cvp call(default: 1)')
    parser.add_argument("--optimize-strategy", type=bool, default=False,
                        help="Optimize Colattice Strategy or not")
    
    args = parser.parse_args()

    colattice_test(
        min_dim=args.min_dim,
        max_dim=args.max_dim,
        step=args.step,
        tours=args.tours,
        approx_factor=args.approx_factor,
        threads=args.threads,
        default_sieve=args.default_sieve,
        max_sample_times=args.max_sample_times,
        batch_size = args.batch_size,
        len_bound = args.len_bound,
        optimize_strategy = args.optimize_strategy
    )
