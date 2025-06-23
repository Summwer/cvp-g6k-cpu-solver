from math import exp, floor, pi, e, log2, ceil, log, lgamma, sqrt
from math import log as ln
from mpmath import mp
from scipy.special import ndtr
# from gate_count.eprint_2019_1161.cost import popcount_grover_iteration_costf

def compute_delta(k):
    """Computes delta from the block size k. Interpolation from the following
    data table:
    Source : https://bitbucket.org/malb/lwe-estimator/
    src/9302d4204b4f4f8ceec521231c4ca62027596337/estima
    tor.py?at=master&fileviewer=file-view-default
    :k: integer
    estimator.py table:
    """

    small = {0: 1e20, 1: 1e20, 2: 1.021900, 3: 1.020807, 4: 1.019713, 5: 1.018620,
             6: 1.018128, 7: 1.017636, 8: 1.017144, 9: 1.016652, 10: 1.016160,
             11: 1.015898, 12: 1.015636, 13: 1.015374, 14: 1.015112, 15: 1.014850,
             16: 1.014720, 17: 1.014590, 18: 1.014460, 19: 1.014330, 20: 1.014200,
             21: 1.014044, 22: 1.013888, 23: 1.013732, 24: 1.013576, 25: 1.013420,
             26: 1.013383, 27: 1.013347, 28: 1.013310, 29: 1.013253, 30: 1.013197,
             31: 1.013140, 32: 1.013084, 33: 1.013027, 34: 1.012970, 35: 1.012914,
             36: 1.012857, 37: 1.012801, 38: 1.012744, 39: 1.012687, 40: 1.012631,
             41: 1.012574, 42: 1.012518, 43: 1.012461, 44: 1.012404, 45: 1.012348,
             46: 1.012291, 47: 1.012235, 48: 1.012178, 49: 1.012121, 50: 1.012065}

    if k != round(k):
        x = k - floor(k)
        d1 = compute_delta(floor(k))
        d2 = compute_delta(floor(k) + 1)
        return x * d2 + (1 - x) * d1

    k = int(k)
    if k < 50:
        return small[k]
    else:
        delta = GH_sv_factor_squared(k)**(1. / (2. * k - 2.))
        return delta


def GH_sv_factor_squared(k):
    return ((pi * k)**(1. / k) * k / (2. * pi * e))



def bkzgsa_gso_len(logvol, i, d, beta=None, delta=None):
    '''
    Return log(GS-lengths) based on GSA. 
    '''
    if delta is None:
        delta = compute_delta(beta)

    return delta**(d - 1 - 2 * i) * exp(logvol / d)
    

#theo d4f 2
def dims4free(beta):
    return floor(beta * ln(4./3.) / ln(beta/(2*pi*exp(1))))


def gaussian_heuristic(log_rr):
    """
    Return norm of shortest vector as predicted by the Gaussian heuristic.

    :param log_rr: vector of ln(norm) Gram-Schmidt norms

    """
    n = len(list(log_rr))
    log_vol = sum(log_rr)
    log_gh =  1./n * (log_vol - ball_log_vol(n))
    return exp(log_gh)

def ball_log_vol(n):
    """
    Return volume of `n`-dimensional unit ball

    :param n: dimension

    """
    return (n/2.) * log(pi) - lgamma(n/2. + 1)



agps20_gate_data = {
          64  :42.5948446291284,   72  :44.8735917172503, 80  :47.4653141889341, 88  :50.0329479433691, 96  :52.5817667347844,  
          104 :55.1130237325179,   112 :57.6295421450947, 120 :60.133284108578,  128 :62.1470129451821, 136 :65.4744488064273, 
          144 :67.951405476229,    152 :70.0494944191399, 160 :72.50927387359,   168 :74.9619105412039, 176 :77.4100782579645, 
          184 :79.3495443657483,   192 :81.7856479853679, 200 :84.2178462414349, 208 :86.646452845262,  216 :89.0717383389617, 
          224 :91.4939375786565,   232 :93.9132560751063, 240 :96.3298751307529, 248 :98.7439563146036, 256 :101.155644837658, 
          264 :104.091650357302,   272 :106.500713866161, 280 :108.907671199501, 288 :111.312627066864, 296 :113.715679081585, 
          304 :116.11691871212,    312 :118.516432037545, 320 :120.914300351043, 328 :123.310600632063, 336 :125.705405925853, 
          344 :128.098785623819,   352 :130.490805751072, 360 :132.881529104042, 368 :135.271015458153, 376 :137.659321707881, 
          384 :140.046501985502,   392 :142.432607773479, 400 :144.817688009257, 408 :147.201789183958, 416 :149.584955436701, 
          424 :151.967228645918,   432 :154.348648518547, 440 :156.729252677678, 448 :159.109076748918, 456 :161.488154445581, 
          464 :163.866517652676,   472 :166.24419650959,  480 :168.621219491327, 488 :170.997613488119, 496 :173.373403883249, 
          504 :175.748614628914,   512 :178.123268319974, 520 :180.931640474467, 528 :183.305745118107, 536 :185.679338509895, 
          544 :188.052439374005,   552 :190.425065356218, 560 :192.797233085084, 568 :195.168958230518, 576 :197.540255559816, 
          584 :199.911138991095,   592 :202.281621644196, 600 :204.651715889082, 608 :207.02143339179,  616 :209.390785157985, 
          624 :211.759781574203,   632 :214.128432446848, 640 :216.496747039019, 648 :218.864734105257, 656 :221.232401924303, 
          664 :223.599758329925,   672 :225.96681073994,  680 :228.333566183483, 688 :230.700031326626, 696 :233.066212496418, 
          704 :235.43211570344,    712 :237.797746662944, 720 :240.163110814653, 728 :242.528213341298, 736 :244.893059185964, 
          744 :247.25765306831,    752 :249.621999499728, 760 :251.986102797502, 768 :254.349967098032, 776 :256.71359636917,  
          784 :259.076994421734,   792 :261.440164920231, 800 :263.803111392861, 808 :266.165837240825, 816 :268.528345816343, 
          824 :270.890640143248,   832 :273.252723321704, 840 :275.614598434176, 848 :277.976268306208, 856 :280.337735739304, 
          864 :282.699003457275,   872 :285.060074111424, 880 :287.420950285349, 888 :289.781634499399, 896 :292.142129214795, 
          904 :294.502436837451,   912 :296.862559721505, 920 :299.222500172584, 928 :301.582260450819, 936 :303.941842773632, 
          944 :306.301249318305,   952 :308.660482224348, 960 :311.019543595679, 968 :313.378435502636, 976 :315.737159983825, 
          984 :318.095719047813,   992 :320.454114674691, 1000:322.8123488175,   1008:325.170423403542,1016:327.52834033558, 
          1024:329.886101492934
          }


# Return log2 of the number of gates for FindAllPairs according to AGPS20
def agps20_gates(beta_prime):
    k = beta_prime / 8
    if k != round(k):
        x = k - floor(k)
        d1 = agps20_gates(8*floor(k))
        d2 = agps20_gates(8*(floor(k) + 1))
        return x * d2 + (1 - x) * d1
    return agps20_gate_data[beta_prime]

# Return log2 of the number of vectors for sieving according to AGPS20
def agps20_vectors(beta_prime):
    k = round(beta_prime)
    N = 1./caps_vol(beta_prime, pi/3.)
    return log2(N)



# Function C from AGPS20 source code
def caps_vol(d, theta, integrate=False, prec=None):
    """
    The probability that some v from the sphere has angle at most theta with some fixed u.

    :param d: We consider spheres of dimension `d-1`
    :param theta: angle in radians
    :param: compute via explicit integration
    :param: precision to use

    EXAMPLE::

        sage: C(80, pi/3)
        mpf('1.0042233739846629e-6')

    """
    prec = prec if prec else mp.prec
    with mp.workprec(prec):
        theta = mp.mpf(theta)
        d = mp.mpf(d)
        if integrate:
            r = (
                1
                / mp.sqrt(mp.pi)
                * mp.gamma(d / 2)
                / mp.gamma((d - 1) / 2)
                * mp.quad(lambda x: mp.sin(x) ** (d - 2), (0, theta), error=True)[0]
            )
        else:
            r = mp.betainc((d - 1) / 2, 1 / 2.0, x2=mp.sin(theta) ** 2, regularized=True) / 2
        return r



list_decoding = {"classical_matzov22": (0.29613500308205365,20.387885985467914), "classical_agps20": (0.2988026130564745, 26.011121212891872), "list_decoding-g": (0.26600114174341505, 23.440974518186337), "list_decoding-naive_quantum": ( 0.2632557273632713, 15.685687713591548) }



def popcount_costf(d, metric = "naive_classical"):
    if(metric ==  "naive_classical"):
        return 1
    elif(metric == "classical_matzov22" or metric == "classical_agps20"):
        n = 1
        while n < d:
            n = 2 * n
            
        ell = mp.ceil(mp.log(n, 2))
        gates = n + (n - ell - 1)*5 + ell
        return gates

    #else:
    #    return popcount_grover_iteration_costf()["gates"]
    

def svp_costf(n, metric = "naive_classical"):
    """
    Return time cost for one svp call.
    """
    if(metric ==  "naive_classical"):
        return 0.292*n
    elif(metric == "classical_matzov22"):
        return list_decoding["classical_matzov22"][0] * n +   list_decoding["classical_matzov22"][1]
    elif(metric == "classical_agps20"):
        return list_decoding["classical_agps20"][0] * n +   list_decoding["classical_agps20"][1]
    elif(metric == "g"):
        return list_decoding["list_decoding-g"][0] * n +   list_decoding["list_decoding-g"][1]
    elif(metric =="naive_quantum"):
        return list_decoding["list_decoding-naive_quantum"][0] * n +   list_decoding["list_decoding-naive_quantum"][1]
        

def bkz_cost(d, beta, metric = "naive_classical", consider_d4f = False):
    """
    Return time cost for one bkz call
    """
    if(consider_d4f):
        f = dims4free(beta)
        beta_prime = floor(beta - f)
    else:
        beta_prime = beta
   
    return svp_costf(beta_prime, metric = metric) + log2(d-beta)
   

C = 1./(1.- 2**(-.292))
def probkz_cost(d, beta, metric = "classical_agps20", consider_d4f = False):
    return log2(C)+ bkz_cost(d, beta, metric = metric, consider_d4f = consider_d4f)









norm_prob_table = {
    0.00: 0.000000, 0.01: 0.007980, 0.02: 0.015960, 0.03: 0.023939, 0.04: 0.031918,
    0.05: 0.039896, 0.06: 0.047873, 0.07: 0.055849, 0.08: 0.063824, 0.09: 0.071798,
    0.10: 0.079770, 0.11: 0.087741, 0.12: 0.095710, 0.13: 0.103677, 0.14: 0.111642,
    0.15: 0.119605, 0.16: 0.127566, 0.17: 0.135524, 0.18: 0.143480, 0.19: 0.151433,
    0.20: 0.158519, 0.25: 0.197413, 0.30: 0.235823, 0.40: 0.310843, 0.50: 0.382925,
    0.60: 0.451494, 0.70: 0.516073, 0.80: 0.576289, 0.90: 0.631880, 1.00: 0.682689,
    1.10: 0.728668, 1.20: 0.769861, 1.30: 0.806399, 1.40: 0.838487, 1.50: 0.866386,
    1.60: 0.890401, 1.70: 0.910869, 1.80: 0.928139, 1.90: 0.942567, 1.96: 0.950004,
    2.00: 0.954500, 2.10: 0.964272, 2.20: 0.972193, 2.30: 0.978552, 2.40: 0.983605,
    2.50: 0.987581, 2.58: 0.990120, 2.60: 0.990677, 2.70: 0.993066, 2.80: 0.994890,
    2.90: 0.996273, 3.00: 0.997300, 3.10: 0.998048, 3.20: 0.998577, 3.29: 0.999000,
    3.30: 0.999033, 3.40: 0.999337, 3.50: 0.999535, 3.60: 0.999664, 3.70: 0.999751,
    3.80: 0.999809, 3.90: 0.999851, 3.99: 0.999936
}




def inf_norm_prob(dim, inf_bound, distance):
    """
    dim: lattice dimension
    inf_bound: infinity norm bound
    distance: distance ln(||w-t||) in Eulidean norm
    
    return log2(prob)
    """

    # T = norm(loc=0, scale=1) 
    t = sqrt(dim)*inf_bound/ exp(distance)
    t = round(t,1)
    if(t > 3.99):
        prob = 1.
    else:
        prob = norm_prob_table[t] #ndtr(t) - ndtr(-t)  # P(-t ≤ X ≤ t)
    try:
        log_prob = log2(prob) * dim 
    except ValueError:
        return -float("inf")
    
    return log_prob




