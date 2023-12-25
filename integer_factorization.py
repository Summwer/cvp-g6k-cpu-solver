from Large_Integer_Factoring_Complete_Process import generate_fac_relations
from N_factor import factor_N
from math import exp


#1. Generate more than len(P) fac-relations
#N = 768055983587580953(60bits) ,n=40 test

N = 768055983587580953
total_trials = 100000 #总实验个数的上限
bits = 60 #N的bit数
n = 50 #W的维数
l_P = 300 #设置list P长度的上限,防止无法跳出循环，l_P的值应该小于等于p_upper的最大素数的位置。
z = round(exp(n/2)) #第n+1行与前n行的权重系数
C = z*2**(100+bits) #C为要放大的精度


sigma = 2 #高斯采样标准差
#1. generate fac-relations
N,uv_list,P = generate_fac_relations(bits,n,l_P,C,z,sigma,total_trials,N)

#2. factor the N
factor_N(N,uv_list,P)



