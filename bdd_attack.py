from g6k.siever_params import SieverParams
from g6k.siever import Siever
from fpylll import IntegerMatrix, LLL
from timu import timu

# 构造反循环矩阵（negacyclic）
def create_negacyclic_matrix(coeffs, n):
    """
    构造一个 n×n 的 negacyclic matrix：模 x^n + 1 的卷积矩阵
    """
    F = IntegerMatrix(n, n)
    for i in range(n):
        for j in range(n):
            idx = (i - j) % n
            sign = -1 if (i - j) < 0 else 1
            F[i, j] = sign * coeffs[idx]
    return F




def lattice_basis(A, q, m=None):
    """
    Construct primal lattice basis for LWE challenge
    ``(A,c)`` defined modulo ``q``.

    :param A: LWE matrix
    :param c: LWE vector
    :param q: integer modulus
    :param m: number of samples to use (``None`` means all)

    """
    if m is None:
        m = A.nrows
    elif m > A.nrows:
        raise ValueError("Only m=%d samples available." % A.nrows)
    n = A.ncols
    B = IntegerMatrix(m+n, m+n)
    for i in range(m):
        for j in range(n):
            B[j, i] = A[i, j]
            B[j, j+n] = 1
        B[i+n, i] = q
        # B[-1, i] = c[i]
    # B[-1, -1] = 1
    
    B = LLL.reduction(B)
    #assert(B[:n] == IntegerMatrix(n, m))
    #B = B[n:]
    print(B)
    return B



def RMLWE_bdd_attack(saitinumber,  LWE_instance_type ="RLWE", params = None):
    params = SieverParams(gpus=2, threads = 32)
    
    print("第%d题" %saitinumber)
    
    n, q, m, a, t, target_norm = timu(saitinumber)
    A = create_negacyclic_matrix(a, n)
    B = lattice_basis(A,q)
    #print(B)
    #c = t[:m]
    #g6k = Siever(B, params)
    
if __name__ == '__main__':
    for saitinumber in [1]:
        RMLWE_bdd_attack(saitinumber)
