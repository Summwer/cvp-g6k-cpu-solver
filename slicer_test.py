from g6k.siever_params import SieverParams
from math import sqrt
from g6k.utils.util import load_cvp_instance
from randslicer_with_d4f import randslicer_with_d4f
import argparse
from g6k.utils.cli import parse_args


def slicer_test(min_dim, max_dim, step, tours, approx_factor, 
               threads, default_sieve, max_sample_times, batch_size, betamax):
    """Run CVP solver tests with given parameters"""
    print("{0: <10} | {1: <10} | {2: <15} | {3: <15} | {4: <20} | {5: <15} | {6: <15} | {7: <20} | {8: <15} | {9: <15} | {10: <15}".format(
        "dim", "index", "batch size", "sample times", "max sample times", 
        "T_sieve (sec)", "T_slicer (sec)", "len_bound", "required len_bound", "db_size", "satisfied vectors"))
    
    for n in range(min_dim, max_dim, step):
        for index in range(1, tours+1):
            
            
            
            A, ts= load_cvp_instance(n, batch_size = batch_size, betamax=betamax)
            
            #sp = SieverParams()
            
            params = SieverParams(threads=threads, 
                                default_sieve=default_sieve, reserved_db_size= max((2 * sqrt(4/3.) ** n), 2000))#, saturation_ratio = 2.)
            
            max_sample_times = max_sample_times
            _, _, ws, _, sample_times, T_pump, T_slicer, db_size, gh = randslicer_with_d4f(
                A, ts, params, len_bound=approx_factor**2, max_sample_times=max_sample_times,kappa=0, blocksize = n, f=0)
            
            dt_min = None
            dt_max = None
            for i in range(len(ts)):
                dt = sqrt(sum([(ws[i][j] - ts[i][j])**2 for j in range(len(ts[0]))]))
                if(dt_min is None or dt < dt_min):
                    dt_min = dt 
                if(dt_max is None or dt > dt_max):
                    dt_max = dt 
            
            if(len(ts) == 1):
                actual_len_bound = round(dt_min/gh,3)
            else:
                actual_len_bound = str( round(dt_min/gh,3))+"~" + str( round(dt_max/gh,3)) #str(round(dt_min,3))+"~"+str(round(dt_max,3))
            print("{0: <10} | {1: <10} | {2: <15} | {3: <15} | {4: <20} | {5: <15} | {6: <15} | {7: <20} | {8: <15} | {9: <15} | {10: <15}".format(
                    n, index, batch_size, sum(sample_times), max_sample_times, 
                    round(T_pump,4), round(T_slicer,4), 
                    actual_len_bound, round(approx_factor,3),
                    db_size, int(.5 * params.saturation_ratio * params.db_size_base ** n)))
            

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
        betamax = args.betamax
    )
