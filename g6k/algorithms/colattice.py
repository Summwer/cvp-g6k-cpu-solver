from .cvpump import cvpump
from g6k.utils.stats import dummy_tracer
from copy import deepcopy
import numpy as np
from fpylll import IntegerMatrix


#We should find the close vector from pi_ln(L) to pi_l1(L) where l1<...<ln
def colattice(g6k, t, blocksizes, sieve_start_index = None, xs = None, yl = None, verbose = True):
    """
    Aim to find the closest vector to target vector t.
    
    
    
    Otherwise, sieve_start_index should be d.
    
    t: target vector 
    
    
    Case: If we want to find the closest vector to a projected vector, the coefficients of [sieve_start_index, d] is fixed. Then we should set the following parameters: 
    sieve_start_index: should be < d.
    xs: the current coefficients w.r.t lattice basis B.
    yl: the current coefficients w.r.t gso basis B*. 
    """
    # if(ee is None):
    
    
    d = g6k.full_n
    # if(pv is None):
    #     pv = [0]*d
    
    
    pt = [0]*d
    pw = [0]*d
    w = [0]*d
    if sieve_start_index is None:
        sieve_start_index = d
    llb = sieve_start_index
    
    if(xs is None):
        xs = [0]*d
        ee = deepcopy(t)
    else:
        mat_x = IntegerMatrix(1,g6k.full_n)
        for i in range(g6k.full_n):
            mat_x[0,i] = int(xs[i])
        # print("xs: ", xs)
        res =  mat_x * g6k.M.B
        print(res)
        ee = tuple(res[0])
        # print("ee:",ee)
    # g6k.initialize_local(0,sieve_start_index,g6k.full_n)
    
    # print([g6k.M.get_r(i,i) for i in range(g6k.full_n)])
    # print(g6k.ll, g6k.l, g6k.r)
    # # print("res: ", res)
    # yl = g6k.M.from_canonical(tuple(res))
    # print("yl: ", yl)

    # print("input xs:", xs)
    # if(sieve_start_index<d):
    #     g6k.initialize_full_cv(v=xs, yl = yl)
    # else:
    #     g6k.initialize_full_cv()
    
    
    # print("input coeff(ee):", [round(_,2) for _ in np.linalg.solve(np.array(list(g6k.M.B)).transpose(), ee)])
        
    # print("d: ", d)
    
    # plus_pv = False
    
    
    # print("ee:", ee)
    
    # print("t:", t)
    
    enter_babai = False
    for i in range(len(blocksizes)-1,-1,-1):
        blocksize = blocksizes[i]
        llb -= blocksize
        if(verbose):
            print("Find close vector on [%d,%d]" %(llb, llb+blocksize))
        if blocksize == 1:
            # if(not plus_pv):
            #     t = [t[i]+ pv[i] for i in range(d)]
            #     plus_pv = True
            
            if(not enter_babai or i==len(blocksizes)-1):
                g6k.initialize_local(0,llb+1, d)
                t_yr =  [round(_,2) for _ in g6k.M.from_canonical(tuple(ee))]
                for i in range(d):
                    if(i<=llb+1):
                        t_yr[i] = 0
                    else:
                        break
                t = g6k.M.to_canonical(tuple(t_yr))
                print("input t:",[round(_,2) for _ in g6k.M.from_canonical(tuple(t))])
                
                
                g6k.initialize_full_cv(v=xs, yl = t_yr)
                
                # while(llb+1<d-50):
                #     g6k.extend_left(1)
                # g6k.initialize_local(0, llb, d)
                # g6k.initialize_target_vector(t,ll = llb+1)
                g6k.initialize_target_vector(t)#,ll = llb)
                enter_babai = True

            
                # print("xs:", xs)
                
                
               
                
            
            if(verbose):
                print("Start cvp_extend_left.")
            
            g6k.cvp_extend_left(1)
            # print(i,sum(blocksizes), d, g6k.ll,g6k.l,g6k.r)
            # print("t(in babai):",t)
            # print("coeff_t(in babai):")
            # print([round(_,2) for _ in g6k.M.from_canonical(tuple(t))]) 
            pw, w, x = g6k.get_cv()
            
            
            # full_pw, w, full_x = g6k.get_full_cv()
            
            
            # print("ee:",ee)
            
            
            # print(g6k.ll, g6k.l, g6k.r)
            # t = pw 
            
            # ee = t
            # ee = [ee[i] - w[i] for i in range(g6k.full_n)]
            # ee = tuple([t[i] - w[i] for i in range(d)])
            # print("coeff(w) in babai: ", [round(_,2) for _ in np.linalg.solve(np.array(list(g6k.M.B)).transpose(), w)]) 
            # print("coeff - ee:", [round(_,2) for _ in np.linalg.solve(np.array(list(g6k.M.B)).transpose(), ee)])
            # print("coeff(ee) (in babai):", [round(_,2) for _ in np.linalg.solve(np.array(list(g6k.M.B)).transpose(), ee)])
            
        
            
            # print("full_pw:", full_pw)
            
            # print("full_x: ", full_x)
            
            
            # print("coeff(t):", [round(_,2) for _ in g6k.M.from_canonical(tuple(t))])
            
            
            # print("coeff(ee):", [round(_,2) for _ in g6k.M.from_canonical(tuple(ee))])
            
            
            
            xs[g6k.l] = x[g6k.l]
            # print("dim=1 xs:", xs)
            
            # print("xs:", xs)
            # print("x:", x)
            
            
            # mat_x = IntegerMatrix(1,g6k.full_n)
            # for i in range(g6k.full_n):
            #     mat_x[0,i] = int(xs[i])
            # # print("xs: ", xs)
            # res =  mat_x * g6k.M.B
            # # print(res)
            # ee = tuple(res[0])
            
            
            ee = w 
            
            # print("coeff(ee) from B*:", [round(_,2) for _ in g6k.M.from_canonical(tuple(ee))])
            
            # raise ""
            
            # raise ""
        else:
            # if(plus_pv):
            #     t = [t[i] - pv[i] for i in range(d)]
            #     plus_pv = False
                
            f = 0
            if(verbose):
                #max(blocksize,50)
                print("Start cvpump-{%d,%d,%d}" %(max(llb,0),blocksize,f))
            t = tuple([t[i] - pt[i] + pw[i] - w[i] for i in range(d)])
            # t = tuple([t[i] - w[i] for i in range(d)]) 
            # print("t:", t)
            pt, pw, w, x = cvpump(g6k,t,dummy_tracer,max(llb,0),blocksize, f,verbose=verbose)
            
            # print("pt:", [ round(_,2) for _ in pt])
            print("projected norm in slicer:", sum([(pt[i]-pw[i])**2 for i in range(len(pw))]))
            # print("t from canonical:" )
            
            # # print("t after cvpump:", pt)
            # print("pw:")
            # print([round(_,2) for _ in g6k.M.from_canonical(tuple(pw))]) 
            # print("coeff(pt):")
            # print([round(_,2) for _ in g6k.M.from_canonical(tuple(pt))]) 
            ee = tuple([ee[i] - w[i] for i in range(d)]) #Final ee = t - sum(w)
            
            # print("coeff(ee) from B:", [round(_,2) for _ in np.linalg.solve(np.array(list(g6k.M.B)).transpose(), ee)])
            
            
            # print("coeff(ee) from B*:", [round(_,2) for _ in g6k.M.from_canonical(tuple(ee))])
            
            
            
            # print("coeff(w):", [round(_,2) for _ in np.linalg.solve(np.array(list(g6k.M.B)).transpose(), w)])
            
            
            
            
            # print("coeff(ee) after slicer:", [round(_,2) for _ in g6k.M.from_canonical(tuple(ee))])
            
            # print("coeff_t:", [round(_,2) for _ in g6k.M.from_canonical(tuple(t))])
            
            
            
        
            
            # print("ee: ", ee)
            xs =[xs[i] + x[i] for i in range(d)]
            enter_babai = False
            
            
            
            # pw, w, x = g6k.get_cv()
            
            # print("-- \nx: ", x)
            
            
            # print("dim>1 xs:", xs)
    # print("norm_ee:", sum([_**2 for _ in ee]))
    
    
    print(ee)
    
    
    mat_x = IntegerMatrix(1,g6k.full_n)
    for i in range(g6k.full_n):
        mat_x[0,i] = int(xs[i])
    res =  mat_x * g6k.M.B
    print(res)
    print(sum([res[0,i]**2  for i in range(g6k.full_n)]))
    print(sum([_**2 for _ in ee]))
         
    return ee,xs

