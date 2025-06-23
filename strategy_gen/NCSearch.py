from math import exp, log,floor, log2, sqrt
from .util import dims4free, gaussian_heuristic, inf_norm_prob, svp_costf, popcount_costf
#from cost import popcount_grover_iteration_costf, compose_k_sequential


def construct_avg_strategy( d, k ):
    Psi = [d // k] * k
    # Adjust for remainder
    remainder = d % k
    if remainder > 0:
        for i in range(remainder):
            Psi[i] += 1
    return Psi


class NC_strategy:
    def __init__(self, min_block = 50, metric = "naive_classical", consider_d4f = False, optimize_strategy = False, batch_size = 0):  
        #For parameter setting
        self.min_block = min_block #the minimal block size for each strategy
        self.metric = metric 
        self.consider_d4f = consider_d4f
        self.optimize_strategy = optimize_strategy #Choose original/optimized colattice strategy
        self.inf_bound = None
        self.dist_bound = None
        self.batch_size = batch_size #log(batch_size) for bach-CVP

    def nc_est(self, log_rr, dist_bound = None, inf_bound = None):
        """
        log_rr: log(gs-lengths)
        distance_bound: log(Euclidean distance bound)
        inf_bound: Infinity distance bound
        """
        if(dist_bound is not None):
            self.dist_bound = dist_bound
            self.gamma1 = sqrt(4/3.)
            self.gamma2 = sqrt(4/3.)
            
            
        elif(inf_bound is not None):
            self.dim = len(log_rr)
            self.inf_bound = inf_bound
            self.gamma1 = sqrt(4/3.)
            self.gamma2 = sqrt(4/3.)
            
            
        if(self.optimize_strategy ):
            return self.NC_strategy_explorer(log_rr)
        else:
            return self.Original_NC_strategy(log_rr)
        
    
    
    def compute_dist(self, log_rr, Psi):
        """
        gamma: the approximate factor for beta1 strategy 
        return log(dist)
        """
        
        
    
        
        log_Delta = 0
        l = 0
        
        for beta in Psi:
            r = l + beta
            if(l == 0):
                log_Delta = log (self.gamma1 * gaussian_heuristic(log_rr[l:r]))
            else:
                log_Delta = log (exp(log_Delta*2.) + (self.gamma2 **2) * gaussian_heuristic(log_rr[l:r])**2)/ 2.
            l += beta
                
        return log_Delta


        

    def cvp_cost(self,beta, first_beta = False):
        #alpha = sqrt(4/3.), if gamma1 or gamma2 = sqrt(4/3.), then prob = 0; if gamma1 or gamma2 = 1. , then prob = -0.149780140929454
        
        
        prob = 0. #-1.6017132519074588e-16 * beta
        if(first_beta):
            if(self.inf_bound is not None):
                #Then we should consider the success probability of infinity norm
                prob = inf_norm_prob(self.dim, self.inf_bound, self.distance) 
            elif(self.gamma1 == 1.):
                prob = -0.149780140929454 * beta
        else:
            if(self.gamma2 == 1.):
                prob = -0.149780140929454 * beta
                
            
            
        if(self.consider_d4f):
            f = dims4free(beta)
            beta_prime = floor(beta - f)
        else:
            beta_prime = beta
            
        try:
            if self.metric in ["naive_classical", "classical_agps20", "classical_matzov22"]:
                return log2(pow(2,svp_costf(beta_prime, metric = self.metric) - self.batch_size ) +  (pow(2,0.0850827230910432 * beta_prime - prob + log2(popcount_costf(beta_prime, metric = self.metric)))))
            # else:
                
            #     n = 1
            #     while n < beta_prime:
            #         n = 2 * n

            #     k = int(1/3. * (n - 1))
            #     N = (4/3.) ** (beta_prime/2.)
            #     cvp_look_cost = popcount_grover_iteration_costf(N, n, k, self.metric)
            #     cvp_looks_per_bucket = N ** (1 / 2.0)
            #     cvp_search_one_cost = compose_k_sequential(cvp_look_cost, cvp_looks_per_bucket).gates
                
                
            #     return log2(pow(2,svp_costf(beta_prime, metric = self.metric) - self.batch_size) +  pow(2,log2(cvp_search_one_cost) - prob))
                
        except OverflowError:
            return float("inf")


    def colattice_cost(self,Psi):
        T = None
        for beta in Psi:
            if(T is None):
                T = self.cvp_cost(beta, first_beta=True)
            else:
                T = log2(2**T+ 2**self.cvp_cost(beta))
                
        return T


    def strategy_search_loop(self, Psi_min, T_min, fixed_blocks, k, beta_max, beta_min, log_rr):
        """
        d: Here d = dim - sum(fixed_blocks)
        """
        dim = len(log_rr)
        upper_bound = beta_max 
        lower_bound = beta_min 
        
        
        
        while(upper_bound >= lower_bound):
            # print("[ %d, %d ]" %(lower_bound, upper_bound))
            # print("len(fixedblocks) = ", len(fixed_blocks))
            # print("k = ", k)
            i = (upper_bound + lower_bound) >> 1
                
            d = dim - sum(fixed_blocks)
            it = len(fixed_blocks)
            remain_k = k - it
            le_size = (d-remain_k*self.min_block)//(i-self.min_block)
            
            if(remain_k == 1 or le_size <= 0 or d // k <= 0):
                return Psi_min,T_min
            else:
                Current_Psi_min  = fixed_blocks + [i] + construct_avg_strategy(d-i, remain_k-1)

                if((d-remain_k*self.min_block)%(i -self.min_block)):
                    Current_Psi_max  = fixed_blocks + [i]*le_size  + [(d-remain_k*self.min_block)%(i -self.min_block) + self.min_block] + [self.min_block]*(remain_k - le_size - 1)
                else:
                    Current_Psi_max  = fixed_blocks + [i]*le_size + [self.min_block]*(remain_k - le_size)
            
            if(Current_Psi_min == Current_Psi_max):
                return Psi_min,T_min
            
            assert(sum(Current_Psi_min) == dim)
            assert(sum(Current_Psi_max) == dim)
            assert(len(Current_Psi_max) == k)
            assert(len(Current_Psi_min) == k)           
            
            
                        
            current_min_Dist = self.compute_dist(log_rr, Current_Psi_max)
            current_max_Dist = self.compute_dist(log_rr, Current_Psi_min)
                
            if(self.dist_bound is not None):
                if(current_min_Dist > self.dist_bound):
                    lower_bound = i + 1
                    continue
                elif(current_max_Dist<=self.dist_bound):
                    upper_bound = i
                    # Evaluate and compare
                    self.distance = current_max_Dist
                    current_T = self.colattice_cost(Current_Psi_min)
                    if current_T < T_min:
                        Psi_min = Current_Psi_min
                        T_min = current_T
                            #return Psi_min, T_min
                    if(upper_bound == lower_bound):
                        return Psi_min, T_min
                else:
                    Psi_min, T_min = self.strategy_search_loop(Psi_min,T_min, fixed_blocks + [Current_Psi_min[it]], k, Current_Psi_max[it+1], Current_Psi_min[it+1],log_rr)

                    upper_bound = i 
                    
                    if(upper_bound == lower_bound):
                        return Psi_min, T_min
                    
            elif(self.inf_bound is not None):
                self.distance = current_max_Dist
                current_T1 = self.colattice_cost(Current_Psi_min)
                
                self.distance = current_min_Dist
                current_T2 = self.colattice_cost(Current_Psi_max)
                
                if(current_T2 < current_T1):
                    tmp = Current_Psi_min
                    Current_Psi_min = Current_Psi_max
                    Current_Psi_max =  tmp
                    
                    tmp = current_T1
                    current_T1 = current_T2
                    current_T2 = tmp
                # prob = inf_norm_prob(self.dim, self.inf_bound, self.distance) 
                # prob_min = inf_norm_prob(self.dim, self.inf_bound, self.compute_dist(log_rr, Psi_min)) 
                
                
                
                # if(True):
                #     print("----------------(1)--------------------")
                #     print("[ %d, %d ]" %(lower_bound, upper_bound))
                #     print("Tmin = ", T_min)
                #     print("Psi_min = ", Psi_min)
                    
                #     print("Current_Psi_min = ", Current_Psi_min)
                #     print("Current_Psi_max = ", Current_Psi_max)
                #     print("current_T1 = ", current_T1)
                #     print("current_T2 = ", current_T2)
                #     print("prob_min = ", prob_min)
                #     print("P_B = ", prob)
                
                # if(prob < prob_min):
                #     lower_bound = i+1
                if(current_T2 > T_min):
                    lower_bound = i + 1 
                elif(current_T1 < T_min):
                    Psi_min = Current_Psi_min
                    T_min = current_T1
                    # prob_min = prob
                    
                    upper_bound = i 
                elif(current_T2 == current_T1):
                    return Psi_min, T_min
   
                if(current_T1 > T_min or upper_bound == lower_bound):
                    Psi_min, T_min = self.strategy_search_loop(Psi_min,T_min, fixed_blocks + [Current_Psi_min[it]], k, Current_Psi_max[it+1], Current_Psi_min[it+1],log_rr)
                    
                    
                
                
        return Psi_min, T_min


    def NC_strategy_explorer(self, log_rr):
        """
        Cost Estimation for the optimized Nearst-Coltticer Strategy.
        """
        d = len(log_rr)

        Psi_min = [d]
        if(self.inf_bound is not None):
            self.distance = self.compute_dist(log_rr,Psi_min)
        T_min = self.colattice_cost(Psi_min)
        
        beta_max = d
        beta_min = 0
        #Find the size of Psi_min. 
        k0 = 2
        for k in range(2, d + 1):
            # Create Psi with k elements, each initialized to floor(d/k)
            Psi = construct_avg_strategy(d, k )
            
            self.distance = self.compute_dist(log_rr, Psi)
            # print(Psi, self.distance, self.dist_bound)
            if(self.dist_bound is not None and self.distance<self.dist_bound) or self.inf_bound is not None:
                # Evaluate and compare
                current_T = self.colattice_cost(Psi)
                if current_T < T_min:
                    Psi_min = Psi
                    T_min = current_T
                    k0 = k +1
                # print("current_T = ", current_T)
                # print("Psi = ", Psi)
                # prob = inf_norm_prob(self.dim, self.inf_bound, self.distance) 
                # print("prob = ", prob)
                # print("----------------")
            elif(self.dist_bound is not None):
                break    
            if(T_min is not None and current_T == float("inf")):
                break
            
        # k0 = k    
        fixed_blocks = []
        beta_max = Psi_min[0]
        beta_min = construct_avg_strategy(d, k0 )[0]
        # print("beta_max = ", beta_max)
        # print("beta_min = ", beta_min)
        for k in range(k0, d+1):
            if(d-k*self.min_block <= 0):
                break
            Psi_min, T_min = self.strategy_search_loop(Psi_min, T_min, fixed_blocks,  k, beta_max, beta_min, log_rr)
            break
        return Psi_min, T_min


    def Original_NC_strategy(self,log_rr):
        """
        Cost Estimation for the original Nearst-Coltticer Strategy.
        The original Nearst-Coltticer Strategy is: [d//k] * (k-1) + [d%k]
        """
        d = len(log_rr)

        Psi_min = [d]
        if(self.inf_bound is not None):
            self.distance = self.compute_dist(log_rr,Psi_min)
        T_min = self.colattice_cost(Psi_min)
        

        #Find the size of Psi_min. 
        for beta in range(d//2, 5, -1):
            if( (d%beta)> 0 and d%beta <5):
                continue
            
            # Create Psi with k elements, each initialized to floor(d/k)
            if(d%beta):
                Psi = [beta] * (d//beta) + [d%beta]
            else:
                Psi = [beta] * (d//beta)
            
            assert(sum(Psi) == d)
            
            
                
            self.distance = self.compute_dist(log_rr, Psi)
            if(self.dist_bound is not None and self.distance<self.dist_bound) or self.inf_bound is not None:
                # Evaluate and compare
                current_T = self.colattice_cost(Psi)
                if current_T < T_min:
                    Psi_min = Psi
                    T_min = current_T

            elif(self.dist_bound is not None):
                return Psi_min, T_min
            
            if(T_min is not None and current_T == float("inf")):
                break
        return Psi_min, T_min
