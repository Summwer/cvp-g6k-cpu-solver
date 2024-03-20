from .cvpump import cvpump
from g6k.utils.stats import dummy_tracer
from copy import deepcopy


#We should find the close vector from pi_ln(L) to pi_l1(L) where l1<...<ln
def colattice(g6k, t, blocksizes, target_norm=None, verbose = True):
    ee = deepcopy(t)
    d = g6k.full_n
    pt = [0]*d
    pw = [0]*d
    w = [0]*d
    llb = d
    for i in range(len(blocksizes)-1,-1,-1):
        blocksize = blocksizes[i]
        llb -= blocksize
    
        #blocksize = 1, using babai to extend lift. 
        if blocksize == 1 and i == len(blocksizes)-1:
            raise "Haven't implemented this function"
        elif blocksize == 1:
            if(verbose):
                print("Start cvp_extend_left.")
            g6k.cvp_extend_left(1)
            pw, w, x = g6k.get_cv()
            ee = tuple([t[i] - w[i] for i in range(d)])
        else:
            f = 0
            if(verbose):
                #max(blocksize,50)
                print("Start cvpump-{%d,%d,%d}" %(max(llb,0),blocksize,f))
            t= tuple([t[i] - pt[i] + pw[i] - w[i] for i in range(d)])
            pt, pw, w, x = cvpump(g6k,t,dummy_tracer,max(llb,0),blocksize, f,goal_r0=target_norm,verbose=verbose)
            ee = tuple([ee[i] - w[i] for i in range(d)])
            
    return ee