from math import floor,pi,exp
from math import log as ln

def dims4free(beta):
    return floor(beta * ln(4./3.) / ln(beta/(2*pi*exp(1))))



for beta in range(60,90,5):
    print(dims4free(beta))