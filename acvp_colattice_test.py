from g6k.utils.util import load_cvp_instance
from fpylll import *
from g6k.siever import Siever
from cvp_g6k_cpu_solver.g6k.algorithms.colattice_backup import colattice
from math import sqrt


n = 100
A, t = load_cvp_instance(n)
A = LLL.reduction(A)

g6k = Siever(A,None)
coblocksizes = [50,50]
ee, _ = colattice(g6k, t, coblocksizes)

norm_ee = sqrt(sum([_**2 for _ in ee]))
print("ee = ", ee)
print("norm(ee) = ", norm_ee)