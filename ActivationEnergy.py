# Activation Energy

from Cantera import *
from SDToolbox import *
import math
from numpy import *
import scipy as Sci
import scipy.linalg
import numpy as np
import time
from scipy.integrate import odeint

import matplotlib.pyplot as plt

def adiabaticTemp(T, P, mech, q, const):

    gas = importPhase(mech);
    gas.set(T = T, P = P, X = q)
    gas.equilibrate(const)

    return gas.temperature()

def activationEnergy(mech, nH2, P0, Tinit, Tprim):
    ActEne = []
    H2perc = []
    R1 = 8.3144621 # J/mol*K
                     
    for i in range(0,len(nH2)):
        # nH2 = 0.3
        nO2 = (1-nH2[i])*0.21;
        nN2 = (1-nH2[i])*0.79;
        q = 'H2:'+str(nH2[i])+' O2:'+str(nO2)+' N2:'+str(nN2)+'';

        T0 = 0.9*adiabaticTemp(Tinit, P0, mech, q, "HP")

        t0 = cv_CJInd(0, P0, T0, q, mech, 0)
                 
        T1 = T0+Tprim
        
        t1 = cv_CJInd(0, P0, T1, q, mech, 0)
     
        AE = (R1*T0*((-T0*(t1-t0))/(t0*Tprim)+2))/1000

        ActEne.append(AE)
        H2perc.append(nH2[i])

    return ActEne  

'''
#############################################
#mech = 'h2air_highT.cti'
mech = 'h2o2mech_Lietal_2003.cti'
#mech = 'gri30.cti'

#nH2 = [x * 0.01 for x in range(21, 60)]
nH2 = [0.12, 0.13, 0.14, 0.15, 0.2, 0.25, 0.29, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]
#nH2 = [x * 0.01 for x in range(18, 59)]
#R1 = 1.9872041347992353 # cal/mol*K

T1fin = []
P0 = 101325
Tinit = 300
Tprim = 30


AE = activationEnergy(mech, nH2, P0, Tinit, Tprim)

ActEnePfit1 = np.poly1d(np.polyfit(nH2, AE, 6))
ActEnePfit = np.polyfit(nH2, AE,  6)

print AE
print ActEnePfit

x = [0.12, 0.13, 0.14, 0.15, 0.2, 0.25, 0.29, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]
ActEneLine = []
for i in x:
    ActEneLine.append(ActEnePfit[0]*i**6+ActEnePfit[1]*i**5+ActEnePfit[2]*i**4+ActEnePfit[3]*i**3+ActEnePfit[4]*i**2+ActEnePfit[5]*i+ActEnePfit[6])
    #ActEneLine.append(ActEnePfit[0]*i**5+ActEnePfit[1]*i**4+ActEnePfit[2]*i**3+ActEnePfit[3]*i**2+ActEnePfit[4]*i+ActEnePfit[5])
    #ActEneLine.append(ActEnePfit[0]*i**2+ActEnePfit[1]*i**1+ActEnePfit[2])
    

print ActEneLine
plt.plot(nH2, AE,'-', x, ActEneLine, '--', nH2, ActEnePfit1(nH2),'ro')
plt.ylim([0,100])
plt.show()
        
'''
