import math

#Compute gcd(X+1,N) and gcd(X-1,N) to factor N

#Generate a prime list with first n primes in order from small to large.
'''
def first_primes(n):
    p = 1
    P = []
    while len(P) < n:
        p = next_prime(p)
        P += [p]
    return P
'''


#计算X_exponent
def generate_X_exponent(d,U1mU2,n,m):
    X_exponent = [0]*n 
    for j in range(n):
        for i in range(m):
            X_exponent[j]+=d[0][i]*U1mU2[i][j]
    X_exponent = [_/2 for _ in X_exponent]
    return X_exponent


#得到X
def compute_X(X_exponent,d,N,P):
    X = 1
    for i in range(len(P)):
        if(X_exponent[i]>=0):
            X*=P[i]**(X_exponent[i])
        else:
            X*=(inverse_mod(P[i],N))**(-X_exponent[i])
    return X  
    
    
def factor_smooth(x,P):
    fac_coeff = [0]*(len(P)+1)
    if(x!=0):
        y = x
        for i in range(len(P)):
            p = P[i]
            while p.divides(y):
                y /= p
                fac_coeff[i] +=1 
    if(x<0):
        fac_coeff[len(P)]=1 #Denote positive and negtive numbers
    return fac_coeff

#Compute the coefficient matrix of u_i and u_i'
def Compute_U1_U2(uv_list,P):
    U1 = []
    U2 = []
    for i in range(len(uv_list)):
        fac = uv_list[i]
        U1.append(factor_smooth(fac[0], P))
        U2.append(factor_smooth(fac[2], P))
        #print(factor(fac[0]))
        #print(factor_smooth(fac[0], P))
    return U1,U2
    
#Compute the minus result of e_ij-e_ij' for all i and j
def coefficient_minux(U1,U2,P):
    U1mU2 = []
    for i in range(len(U1)):
        row = []
        for j in range(len(U1[i])-1):
            row.append( U1[i][j]-U2[i][j])
        row.append(U2[i][len(P)])
        U1mU2.append(row)
    return U1mU2
 
#Solve the binary vector d to get a square number X
def Solve_D(U1mU2):
    U1mU2_mod_2 = matrix(GF(2),U1mU2).T
    #b = matrix(GF(2),[0]*(U1mU2_mod_2.nrows())).T
    D = U1mU2_mod_2.right_kernel().basis()
    return D

    
    
    
def factor_N(N,uv_list,P):
    l_P = len(P)
    l_uv = len(uv_list)
    #P = first_primes(l_P) #Get list P
    U1,U2 = Compute_U1_U2(uv_list,P) #Calculate Efficients of u_i and u_i' expressed by P.
    U1mU2 = coefficient_minux(U1,U2,P) #e_ij-e_ij'
    D = Solve_D(U1mU2) #solve d to ensure X is square
    #try all solutions in D to find a appropriate X to factor N
    N_factor = []
    for i in range(len(D)):
        d = D[i]
        d = matrix(ZZ,d)
        X_exponent = generate_X_exponent(d,U1mU2,l_P,l_uv)
        X = compute_X(X_exponent,d,N,P)
        if(gcd(X+1,N)!=(1 or N) or gcd(X-1,N)!=(1 or N)):
            N_factor=[gcd(X+1,N),gcd(X-1,N)]
        print(gcd(X+1,N))
        print(gcd(X-1,N))
        break
    return N_factor

