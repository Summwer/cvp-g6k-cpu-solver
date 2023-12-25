from fpylll import *
import sys
from tqdm import tqdm
import time
from random import Random
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from sympy import *
from math import log as ln
from g6k.algorithms.pro_randslicer import pro_randslicer
from g6k.siever import Siever
from g6k.utils.stats import dummy_tracer

'''

My sieving algorithm for sieve fac-relations


Integer factoring --> CVP solver(inf-norm<=lnpn) for the same N 

Matrix = [[f(1)   0    ...   0   ]
          [0     f(2)  ...   0   ]
          [0      0    ...  f(n) ]
          [lnp1 lnp2   ... lnpn ]]
        
target vector = [0,...,0,lnN].T


'''

def svp(B):
    A = IntegerMatrix.from_matrix(B)
    return SVP.shortest_vector(A)


#这个cvp直接就是求的L(B)到k-target vector的距离
def cvp(B,t):
    A = B
    A = LLL.reduction(A)
    t = tuple(t)
    
    g6k = Siever(A)
    print(A.ncols,A.nrows)
    print(g6k.M.B.ncols, g6k.M.B.nrows, g6k.full_n)
    f = 0
    print("GSO precision: ", g6k.M.float_type)
    # t = tuple([t[i] - ]
    close_vector,_ = pro_randslicer(g6k,t,dummy_tracer,f,verbose=True)
    print(close_vector[-1],t[-1],_[-1])
    return close_vector #CVP.closest_vector(A,w)

#Generate a prime list with first n primes in order from small to large.
def first_primes(n):
    p = 1
    P = []
    while len(P) < n:
        p = nextprime(p)
        P += [p]
    return P


#Judge whether could add prime into prime database
def prime_database_addtive_judgement(u,u_,upper,l_P,P,P0,add_size):
    prime_factor_list = sorted(list(set([_[0] for _ in list(factor(u))]+[_[0] for _ in list(factor(u_))])))
    for p in prime_factor_list:
        if(p>upper):
            #print("Add Failure:Since Factor in prime_factor_list is larger than p_upper=%d" %upper)
            return False
    #Determine whether the database will scall up larger than l_P
    Addtive_list = sorted(list(set(prime_factor_list+P)))
    if(len(Addtive_list)>l_P):
    #    print("Add Failure:Since if we add factors into P ,then it will be larger than l_P = %d" %l_P)
        return False
    #When P!=[], primes factoring u_ add into P should be no more than add_size
    Addtive_list_for_u_P = sorted(list(set([_[0] for _ in list(factor(u_))]+P)))
    Addtive_list_for_u_P_0 = sorted(list(set([_[0] for _ in list(factor(u_))]+P0)))
    if(len(Addtive_list_for_u_P_0)-len(P0)>add_size):
        #print("Add Failure:Since addtive number of factors for u_  equal to %d is more than add_size=%d based on both P0"
        #      %(len(Addtive_list_for_u_P_0)-len(P0),add_size))
        return False
    if(P!=[] and len(Addtive_list_for_u_P)-len(P)>add_size):
        #print("Add Failure:Since addtive number of factors for u_  equal to %d is more than add_size=%d based on both P"
        #      %(len(Addtive_list_for_u_P)-len(P0),add_size))
        return False
    return True
        
            
#input u and u',generate a prime database,set an upper factor prime value
def generate_prime_database(u,u_,P):
    prime_factor_list = sorted(list(set([_[0] for _ in list(factor(u))]+[_[0] for _ in list(factor(u_))])))
    Addtive_list = sorted(list(set(prime_factor_list+P)))
    #if(u_<0):
    #    Addtive_list+=[-1]
    #Addtive_list = sorted(Addtive_list)
    return Addtive_list

#test whether the specific number is smooth.
def is_smooth(x, P):
    if(x!=0):
        y = int(x)
        # print(int(y))
        for p in P:
            while(y%p == 0):
                y /= p
        return abs(y) == 1
    if(x==0):
        return False

#calculate the norm of vector
def norm(x,x_p): #x_p is the norm value of x
    norm = 0
    for i in x:
        norm+=abs(i**x_p)
    norm = float(norm**(1/x_p))
    return norm


# Scale up and round 
def sr(x,C,z):
    return round(x*C*z)
    
def Large_Integer_Factoring_reduce_to_kCVP(C,P,n,f,z):
    '''

    Matrix = [           f(i)
          [lnp1 lnp2   ... lnpn ]]

    '''
    # B = IntegerMatrix(n,n+1)
    # for i in range(n):
    #     B[i,i] = f[i]*C
    #     B[i,n] = sr(ln(P[i]),C,z) #跳着选p,每次选第l_P/n*i个p
    # return B
    
    B = IntegerMatrix(n,n)
    for i in range(n):
        B[i,i] = f[i]*C
        B[i,n-1] = sr(ln(P[i]),C,z) #跳着选p,每次选第l_P/n*i个p
    return B



#initialize f， 使生成的随机数服从单边gauss分布，1最多，其它数字逐渐减少。
def generate_f(n,sigma):
    D = DiscreteGaussianDistributionIntegerSampler(sigma) #再sigma上进行高斯采样
    f = []
    #采样时只取>=0的结果
    while(len(f)<n):
        d = D()
        if(d>=0):
            f.append(d+1)
    return f

# Test if a factoring relation was indeed found.
def test_Schnorr(C,N, n,l_P,uv_list,z,sigma,trials, total_trials):
    P0 = first_primes(n)
    P = first_primes(l_P)
    f = generate_f(n,sigma)
    #-----New matrix for k-CVP --------#
    B = Large_Integer_Factoring_reduce_to_kCVP(C,P0,n,f,z) # n*(n+1)
    # target_vector = [0 for _ in range(n)]+[sr(ln(N),C,z)] #dimension is n+1
    target_vector = [0 for _ in range(n-1)]+[sr(ln(N),C,z)] #dimension is n
    #--------------------------------------#
    b = cvp(B,target_vector) #solve the nearest vector of L(M) to target vector.
    e = [b[i]/f[i]/C for i in range(n-1)] #solve coefficient [e1,...en]. 
    e.append(round((b[n-1] - sum([e[i]*sr(ln(P0[i]),C,z) for i in range(n-1)]))/C/ln(P0[n-1])))
    #generate u and v
    u = 1
    v = 1
    if(P==[]):
        for i in range(n):
            # assert e[i] in ZZ
            if e[i] > 0:
                u *= P0[i]**e[i] 
            if e[i] < 0:
                v *= P0[i]**(-e[i])
                
    elif(P!=[]):
        for i in range(n):
            # assert e[i] in ZZ
            if e[i] > 0:
                u *= P0[i]**e[i] 
            if e[i] < 0:
                v *= P0[i]**(-e[i])
    
    u_ = u-v*N
    #Judge_additive
    print(u,u_)
    if((u,v,u_) not in uv_list and is_smooth(u_, P)):
        uv_list.append((u,v,u_))
        print('\r Test trials=%d/%d,len(uv_list)=%d,len(P)= %d，gap= %d. '%(trials,total_trials,len(uv_list),len(P),len(P)-len(uv_list)),end='\n')
    else:
        print('\r Test trials=%d/%d,len(uv_list)=%d,len(P)= %d，gap= %d. '%(trials,total_trials,len(uv_list),len(P),len(P)-len(uv_list)),end='')
        
        
    return is_smooth(u - v*N, P),uv_list

    
def generate_fac_relations(bits,n,l_P,C,z,sigma,total_trials,N):
    start_time = time.time()
    print("[My reduction(with f)]z=%d, C=%d, n=%d,l_P=%d,sigma = %f, %d bits, %d total trials"%(z,C,n, l_P,sigma, bits, total_trials))
    #print('P=',end='')
    #print(first_primes(l_P))
    successes = 0
    #p = random_prime(2**(bits/2), false, 2**(bits/2-1)) #generate a big prime p with length of bits/2
    #q = random_prime(2**(bits/2), false, 2**(bits/2-1)) #generate a big prime q with length of bits/2
    #N = p*q #construct N,fixed N
    #N=499107764672820907790489294293
    #N = 768055983587580953
    print('N=%d.'%N)
    max_C = 0
    max_success = 0
    uv_list = []
    i = 0
    trials = 0
    #We need len(uv\_list)>=len(P)+1
    while( len(uv_list) < l_P+1 and trials < total_trials):
        
        success,uv_list= test_Schnorr(C,N, n,l_P,uv_list,z,sigma,trials,total_trials)  
        successes += success
        
        i = len(uv_list)
        trials +=1
        
        sys.stdout.flush()
        
        pass
    end_time = time.time()
    print("\n %d Factoring Relation found out of %d trials and the valid fac-relations are %d in time %f s"%(successes, trials,len(uv_list),end_time-start_time))
    print('N=%d' %N)
    print('uv_list=',end='')
    print(uv_list)
    print('P=',end='')
    P = first_primes(l_P)
    print(P)
    
    return N,uv_list,P
    
    
    
 

