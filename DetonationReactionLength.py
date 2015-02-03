'''
Calculations of detonation parameters + reaction Lenght
File use zndMZ_reactionLength.py which is based on 
SDToolbox build in Matlab

'''

import numpy as np
import pandas as pd
from SDToolbox import *
from zndMZ_reactionLength import *

# define reaction mechanism
mech = 'h2air_highT.cti'
P1 = 100000;
T1 = 300;

fi = [];
mfH2 = []
mfO2 = []
mfN2 = []

fcsv = open('CJ_state_H2.csv','w')
writeCSV (fcsv, ['T1 [K]', 'P1 [bar]', 'nH2[mol]', 'nO2[mol]', 'nN2[mol]', 'mfH2[-]', 'mfO2[-]', 'mfN2[-]', 'P_cj[bar]', 'T_cj[K]', 'CJ_speed[m/s]', 'ss_af[m/s]', 'ss_ae[m/s]', 'cellsize[mm]'])

# 1 - calculate CJ reaction length, 0 - do not calculate
rCJ = 1
nH2 = [x * 0.01 for x in range(18, 19)]
for i in range(0,len(nH2)):
    nO2 = (1-nH2[i])*0.21;
    nN2 = (1-nH2[i])*0.79;
    q = 'H2:'+str(nH2[i])+' O2:'+str(nO2)+' N2:'+str(nN2)+'';
    fi.append((nH2[i]/(1-nH2[i]))*(4.76/2))
    
    mfH2.append(nH2[i]/(nH2[i]+nO2+nN2))
    mfO2.append(nO2/(nH2[i]+nO2+nN2))
    mfN2.append(nN2/(nH2[i]+nO2+nN2))

    cjSpeed = CJspeed(P1, T1, q, mech, 0);
    [cj_speed,R2] = CJspeed(P1, T1, q, mech, 0);   
    gas = PostShock_eq(cj_speed, P1, T1, q, mech)
    Ps = gas.pressure()/P1
    [ae,af] = equilSoundSpeeds(gas)
    
    # length of the reaction front from the pressure peak to the max Thermicity or max temperature gradient
    deltaInd = zndMZ_reactionLength.znd_noshk(P1, T1, q, mech, rCJ)[1][14]
    # Thermicity
    thermicity = zndMZ_reactionLength.znd_noshk(P1, T1, q, mech, rCJ)[1][7]
    
    writeCSV(fcsv, [T1, P1, nH2[i], nO2, nN2, mfH2[i], mfO2[i], mfN2[i], Ps, gas.temperature(), cj_speed, af, ae, deltaInd[0]])

fcsv.close()   

print "File written to CJ_state_H2.csv"

