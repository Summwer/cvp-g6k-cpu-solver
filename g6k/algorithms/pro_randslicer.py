

def norm(vec):
    return sum([_**2 for _ in vec])

def pro_randslicer(len_bound = 1., max_sample_times = 1000):
  
    _,pro_randslicer.w,_,sample_times  = pro_randslicer.g6k.randslicer(max_sample_times = max_sample_times,len_bound = len_bound)
            

    pro_randslicer.ee = [pro_randslicer.c[i] - pro_randslicer.w[i] for i in range(pro_randslicer.g6k.M.B.ncols)]
    pro_randslicer.norm_ee = norm(pro_randslicer.ee)
            
            
    return pro_randslicer.w, pro_randslicer.ee, sample_times
