from Cantera import *
import sys
import numpy as np
from numpy import *
import math
import matplotlib.pyplot as plt


def cpPolynomials(P1,q,mech,tempRange1,tempRange2,step):

    gas = importPhase(mech);
    # specific heat of fresh mixture
    
    a = tempRange1[0]
    b = tempRange1[1]

    Tv1 = np.arange(a,b,step);
    cp1 = [];
        
    for i in range(a,b,step):
        gas.set(T = i, P = P1, X = q);
        cp1.append(gas.cp_mass());

    CpPoly1 = np.polyfit(Tv1, cp1, 6)

    #print 'Polynomials temperature range:'+str(a)+'-'+str(b)+''
    #print CpPoly1

    
    c = tempRange2[0]
    d = tempRange2[1]

    Tv2 = np.arange(c,d,step);
    cp2 = [];
        
    for i in range(c,d,step):
        gas.set(T = i, P = P1, X = q);
        cp2.append(gas.cp_mass());

    CpPoly2 = np.polyfit(Tv2, cp2, 6)
    
    #print 'Polynomials temperature range:'+str(c)+'-'+str(d)+''
    #print CpPoly2

    print 'Function cpPolynomials teturn list Tv1, cp1, CpPoly1, Tv2, cp2, CpPoly2'
    return [Tv1, cp1, CpPoly1, Tv2, cp2, CpPoly2]

######################################################################
'''
P1 = 100000;
T1 = 298;
q = 'H2:2 O2:1 N2:3.76'
# xCh4 + (1-x)H2 + (1,5x+0,5)O2+3,76
mech = 'h2air_highT.cti'
tempRange1 = [290, 1000]
tempRange2 = [1000, 5000]
step = 10

out = cpPolynomials(P1,q,mech,tempRange1,tempRange2, step)
print out[2]
print out[5]

'''
