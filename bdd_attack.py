from g6k.siever_params import SieverParams
from g6k.siever import Siever
from fpylll import IntegerMatrix, LLL
from timu import timu
from colattice import colattice

def load_inv_mid_matrix(filename):
    """
    Load mimabisai RLWE challenge from file or website.

    :param saitinumber: Number of saiti.

    """
    try:
        data = open(filename, "r").readlines()
    except FileNotFoundError:
        return None
   
    t =  tuple(eval(data[0].replace(" ", ", ")))
    B = eval(",".join([s_.replace(" ", ", ") for s_ in data[1:]]))
    B = IntegerMatrix.from_matrix(B)
    
    return t, B







def store_inv_mid_matrix(saitinumber, B, t_inv):
    #write the mid result of basis       
    filename = 'mimabisai/No-%d-inv-matrix-for-bdd.txt' % (saitinumber)
    fn = open(filename, "w")
    fn.write('[')
    for i in range(B.ncols()):
        if(i< B.ncols()/2):
            fn.write(str(t_inv[i]))
        else:
            fn.write(str(0))
        if i<B.ncols()-1:
            fn.write(' ')
    fn.write(']\n')
    fn.write('[')
    for i in range(B.nrows()):
        fn.write('[')
        for j in range(B.ncols()):
            fn.write(str(B[i][j]))
            if j<B.ncols()-1:
                fn.write(' ')
        if i < B.nrows()-1:
            fn.write(']\n')
    fn.write(']]')
    fn.close()




def RMLWE_bdd_attack(saitinumber, coblocksizes,verbose = False, inv = True):
    params = SieverParams(gpus=2, threads = 32)
    
    
    print("第%d题" %saitinumber)
    
    if(inv):
        t, B = load_inv_mid_matrix('mimabisai/No-%d-inv-matrix-for-bdd.txt' % (saitinumber))
    else:
        t, B = load_inv_mid_matrix('mimabisai/No-%d-matrix-for-bdd.txt' % (saitinumber))
    #n, q, m, a, t, target_norm = timu(saitinumber)
    print(t)
    params = SieverParams(threads = 32)
    g6k = Siever(B, params)
    _ , ee, _ , _, _, _, _= colattice(g6k.M.B, [t], 1., params, blocksizes = coblocksizes, len_bound = 1.,verbose = True)#, target_norm=target_norm)
    ee = ee[0]
    print(ee)
    
if __name__ == '__main__':
    coblocksizes = {1: [32]*4, 5: [64]*4, 7: [64]*8}
    inv = {1: True, 5: False, 7: True}
    for saitinumber in [7]:
        RMLWE_bdd_attack(saitinumber, coblocksizes = coblocksizes[saitinumber],verbose = True, inv= inv[saitinumber])
