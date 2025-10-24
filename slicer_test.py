from g6k.siever_params import SieverParams
from math import sqrt
from g6k.utils.util import load_cvp_instance
from randslicer_with_d4f import randslicer_with_d4f
import argparse
from g6k.utils.cli import parse_args
from strategy_gen.util import dims4free

def slicer_test(min_dim, max_dim, step, tours, approx_factor, 
               threads, default_sieve, max_sample_times, batch_size, betamax, kappa,  consider_d4f, d4f0):
    """Run CVP solver tests with given parameters"""
    print("{0: <10} | {1: <12} | {2: <15} | {3: <15} | {4: <20} | {5: <16} | {6: <15} | {7: <20} | {8: <15} | {9: <18} | {10: <15} | {11: <15}".format(
        "dim", "(ll, l, r)",  "index", "batch size", "sample times", "max sample times", 
        "T_sieve (sec)", "T_slicer (sec)", "len_bound", "required len_bound", "db_size", "satisfied vectors"))
    
    
    
    
    for n in range(min_dim, max_dim, step):
        if(consider_d4f):
            if(d4f0 < 0):
                d4f = dims4free(n-kappa)
            else:
                d4f = d4f0
        else:
            d4f = 0
        
        for index in range(1, tours+1):    
            
            A, ts= load_cvp_instance(n, batch_size = batch_size, betamax=betamax)
            
            #sp = SieverParams()
            
            params = SieverParams(threads=threads, 
                                default_sieve=default_sieve, reserved_db_size= min((2 * sqrt(4/3.) ** (n-kappa-d4f), 2000)))#, saturation_ratio = 2.)
            
            max_sample_times = max_sample_times
            pts, pws, _, _, sample_times, T_pump, T_slicer, db_size, gh = randslicer_with_d4f(
                A, ts, params, len_bound=approx_factor**2, max_sample_times=max_sample_times,kappa=kappa, blocksize = n-kappa, f=d4f)
            
            
            
            dt_min = None
            dt_max = None
            for i in range(len(pts)):
                dt = sqrt(sum([(pws[i][j] - pts[i][j])**2 for j in range(len(pts[0]))]))
                if(dt_min is None or dt < dt_min):
                    dt_min = dt 
                if(dt_max is None or dt > dt_max):
                    dt_max = dt 
            
            if(len(ts) == 1):
                actual_len_bound = round(dt_min/gh,3)
            else:
                actual_len_bound = str( round(dt_min/gh,3))+"~" + str( round(dt_max/gh,3)) #str(round(dt_min,3))+"~"+str(round(dt_max,3))
            
            print("{0: <10} | {1: <12} | {2: <15} | {3: <15} | {4: <20} | {5: <16} | {6: <15} | {7: <20} | {8: <15} | {9: <18} | {10: <15} | {11: <15}".format(
                    n, str((kappa,kappa+d4f,n-kappa)), index, batch_size, sum(sample_times), max_sample_times, 
                    round(T_pump,4), round(T_slicer,4), 
                    actual_len_bound, round(approx_factor,3),
                    db_size, int(.5 * params.saturation_ratio * params.db_size_base ** (n-kappa-d4f))))
            

#Use the simulator like BKZ 2.0
def draw_cvp_bound_simulation():
    return






if __name__ == '__main__':
    description = "Run CVP solver tests with configurable parameters"
    
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
    parser.add_argument('--approx-factor', type=float, default=0.9,
                       help='Approximation factor (default: 0.9)')
    parser.add_argument('--threads', type=int, default=1,
                       help='Number of threads (default: 1)')
    parser.add_argument('--default-sieve', type=str, default="bdgl2",
                       help='Default sieve algorithm (default: bdgl2)')
    parser.add_argument('--max-sample-times', type=int, default=140,
                       help='Maximum sample times (default: 140)')
    parser.add_argument('--batch-size', type=int, default=1,
                       help='number of target vector (default: 1)')
    parser.add_argument('--betamax', type=int, default=55,
                       help='maximal BKZ-beta preprocess')
    parser.add_argument('--kappa', type=int, default=0,
                       help='start index')
    parser.add_argument('--d4f', type=int, default=-1,
                       help='dimension for free value for slicer')
    parser.add_argument('--consider-d4f', type=bool, default=False,
                       help='choose consider dimension for free value for slicer or not')
    
    args = parser.parse_args()

    slicer_test(
        min_dim=args.min_dim,
        max_dim=args.max_dim,
        step=args.step,
        tours=args.tours,
        approx_factor=args.approx_factor,
        threads=args.threads,
        default_sieve=args.default_sieve,
        max_sample_times=args.max_sample_times,
        batch_size = args.batch_size,
        betamax = args.betamax,
        kappa = args.kappa,
        consider_d4f = args.consider_d4f,
        d4f0 = args.d4f
    )
