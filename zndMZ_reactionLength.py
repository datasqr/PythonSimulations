# 1. Add this file to SDToolbox 
# C:\Python27\Lib\site-packages\SDToolbox
# Open _init_.py
# write line: from zndMZ_reactionLength import *
# run the _init_.py (F5)
# enjoy calculations

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


def znd_shk(U1, P1, T1, q, mech, rCJ):
    
    gas1 = importPhase(mech);
    
    gas1.set(T = T1, P = P1, X = q);
    
    gas = PostShock.PostShock_fr(U1, P1, T1, q, mech);
    
    out = znd_detonation(gas,gas1,U1,rCJ);
    
    return out

    
def znd_noshk(P1, T1, q, mech, rCJ):

    plt_num = 1
    cj_speed1 = PostShock.CJspeed(P1, T1, q, mech, plt_num);
    
    cj_speed = cj_speed1[0]
    
    outCJ = znd_shk(cj_speed, P1, T1, q, mech, rCJ)
    
    return [cj_speed, outCJ]


'''
%Set of ODE's to solve ZND Detonation Problem
%
% FUNCTION
% SYNTAX
% dydt = ZNDReactor(t,y,gas,U1,r1,PSC)
% ALWAYS called by and ODE solver, in this case: out = ode15s(@ZNDReactor,tel,y0,options,gas,U1,r1,PSC)
%       tel = time span
%       y0 = initial conditions array
%       options = specified in znd_detonation.m
%       gas = Cantera gas object
%       U1 = Shock Velocity
%       r1 = initial density
%       PSC = Post-shock pressure
%
% INPUT
% t = time
% y = solution array
% gas = Cantera gas object
% U1 = Shock Velocity
% r1 = initial density
% PSC = Post-shock pressure
%
% OUTPUT
% dydt = Array of ODEs to be solved by ode15s
%
% SUBFUNCTION CALLS
% Cantera Functions: set.m, meanMolecularWeight.m, gasconstant.m,
%       density.m, nSpecies.m, netProdRates.m, enthalpies_RT.m,
%       molecularWeights.m, cp_mass.m, soundspeed.m,  
'''

def ZNDReactor(y,t,gas,U1,r1,PSC):

    gas.set(Rho = y[1], Y = y[3:])
    
    wt = gas.meanMolecularWeight();

    rho = gas.density();
    T = (y[0]*PSC/y[1])*(wt/GasConstant);

    gas.set(T = T, Rho = rho, Y = y[3:]); 
    nsp = gas.nSpecies();
            
    # Vectors
    wdot = gas.netProductionRates();

    hs = gas.enthalpies_RT()*GasConstant*T;
    mw = gas.molecularWeights();
        
    # Scalars
    cp = gas.cp_mass();

    c = soundSpeedM(gas);
    #c = equilSoundSpeeds(gas)[0]
    
    U = U1*r1/rho;
    M = U/c;                       #Mach Number
    eta = 1 - M**2;                 #Sonic Parameter

    # % Loop through all of the species,calculate thermicity for this time step
    dykdt = array([0]*nsp, dtype=np.float);
    sum = 0;
    for z in range(0,nsp):
        dykdt[z] = mw[z]*wdot[z]/rho; #Net production rate of species z (kg/sec)
        drdy = -wt/mw[z];             # mean molecular weight / molecular weight of species z
        a = hs[z]/(cp*T*mw[z]);       # factor that accounts for chemical energy change
        sum = sum - (drdy + a)*dykdt[z];
        
    sigmadot = sum;                 #Thermicity

    Pdot = -rho*U**2*sum/eta/PSC;   #Pressure Derivative (Normalized by Post-Shock pressure to reduce size)
    rdot = -rho*sum/eta;           #Density Derivative
    #xdot = U;                      %Distance Derivative
    '''
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        % SET UP COLUMN VECTOR FOR DYDT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        '''
    list3 = array([Pdot, rdot, U]);
    list4 = array([0]*nsp, dtype=np.float)
    dydt = np.concatenate((list3, list4));
        
    # species equations
        
    for i in range(0,nsp):
        dydt[i+3]=mw[i]*wdot[i]/rho;
        
    return dydt

def znd_detonation(gas,gas1,U1, rCJ):

    r1 = gas1.density();
    
    P1 = gas1.pressure();
    T1 = gas1.temperature();
    PSC = gas.pressure();
    rnew = gas.density();
    nsp = gas.nSpecies();

    list1 = array([1, rnew, 0]);
    list2 = array(gas.massFractions());

    y0 = np.concatenate((list1, list2));
    
    tel = np.linspace(0, 0.0001, 10000);

    out = scipy.integrate.odeint(ZNDReactor, y0, tel, args=(gas,U1,r1,PSC),
                                full_output = 0, col_deriv = True,
                                printmessg = 1, rtol = 1.e-5, atol = 1.e-8)

    a,b = out.shape

    i = 1;
    j = 1;
    
    y = array([0]*nsp, dtype=np.float) 

    #Initalize ouput matrices
    outputInd_len_ZND = 0;
    outputInd_time_ZND = 0;
    outputExo_len_ZND = 0;
    outputExo_time_ZND = 0;
    outputTime = array([0]*a, dtype=np.float) ;
    outputDistance = array([0]*a, dtype=np.float) ;
    outputP = array([0]*a, dtype=np.float) ;
    outputT = array([0]*a, dtype=np.float) ;
    outputU = array([0]*a, dtype=np.float) ;
    outputRho = array([0]*a, dtype=np.float) ;
    outputThermicity = array([0]*a, dtype=np.float) ;
    outputM = array([0]*a, dtype=np.float) ;
    outputAf = array([0]*a, dtype=np.float) ;
    outputG = array([0]*a, dtype=np.float) ;
    outputWt = array([0]*a, dtype=np.float) ;
    outputSonic = array([0]*a, dtype=np.float) ;
    outputSpecies = zeros(shape=(a,nsp), dtype=np.float);
    temp_grad = array([0]*a, dtype=np.float);

    # Extract TIME, TEMPERATURE, PRESSURE, and DENSITY arrays from integrator output

    for i in range(0,a):
        
        outputTime[i] = tel[i];

        for m in range(3,nsp+3):
            y[m-3] = abs(out[i,m]);
            
        gas.set(Rho=out[i,1], Y = y);

        wt = gas.meanMolecularWeight();

        Tout = (out[i,0]*PSC/out[i,1])*(wt/GasConstant);
        
        gas.set(T = Tout, P = out[i,0]*PSC, Y = y);
            
        outputP[i] = out[i,0]*PSC/OneAtm;
        outputRho[i] = out[i,1];
        outputDistance[i] = out[i,2];
        outputT[i] = Tout;
        for k in range(1,nsp):
            outputSpecies[i,k] = y[k];


        # Extract WEIGHT, GAMMA, SOUND SPEED, VELOCITY, MACH NUMBER, c^2-U^2,
        # THERMICITY, and TEMPERATURE GRADIENT
 
        den = outputRho[i];
        nsp = gas.nSpecies();

        #Vectors
        wdot = gas.netProductionRates();
        hs = gas.enthalpies_RT()*GasConstant*Tout;
        mw = gas.molecularWeights();

        #Scalars
        cp = gas.cp_mass();
        cv = gas.cv_mass();
        g = cp/cv;
        af = sqrt(g*GasConstant/wt*Tout); #% soundspeed(gas);
        r = gas.density();
            
        U = U1*r1/r;
        M = U/af;                               #Mach Number
        eta = 1 - M**2;                          #Sonic Parameter
        sonic = eta*af**2;

        sum = 0;
        for n in range(1,nsp):
            h = hs[n]/mw[n];
            wd = wdot[n];
            w = mw[n];
            dykdt = w*wd/den;
            drdy = -den*wt/w;
            term = den*h/(cp*Tout);
            sum = sum - (drdy + term)*dykdt;
            
        sigma = sum/den;                        #Thermicity
        Pdot = -U**2*sum/eta/PSC;               #Pressure Derivative
        rdot = -sum/eta;                        #Density Derivative

        # FIND TEMPERATURE GRADIENT
        temp_grad[i] = Tout*(Pdot/out[i,0] - rdot/r);

            
        # Assign ouptput structure
        outputU[i] = U;
        outputThermicity[i] = sigma;
        outputM[i] = M;
        outputAf[i] = af;
        outputG[i] = g;
        outputWt[i] = wt;
        outputSonic[i] = sonic;

    #%Find INDUCTION TIME and LENGTH based on MAXIMUM TEMPERATURE GRADIENT and Maximum
    #%Thermicity
        
    # Procdure to find max value with 
    maxTgrad = max(temp_grad)
    maxThermicity = max(outputThermicity)
    k1 = np.where(temp_grad == max(temp_grad))
    n = np.where(outputThermicity == max(outputThermicity))

    #print 'reaction zone according to MaxGradTemp'
    rZoneGradTemp = outputDistance[k1]
    #print rZoneGradTemp
    
    outputInd_time_ZND = outputTime[n];
    #print 'reaction zone according to MaxThermicity'
    outputInd_len_ZND = outputDistance[n];
    #print outputInd_len_ZND
    
    maxT = max(outputT[:])
    k = np.where(outputT == max(outputT[:]))
    # maxT,k = outputT[:].max(0),outputT[:].argmax(0)

    #print 'reaction zone according to MaxTemp'
    rZoneMaxTemp = outputDistance[k]
    #print rZoneMaxTemp

    
    #############################################################
    #%Find reaction LENGTH based on CJ condition M = 1.0

    if(rCJ == 1):
        Ind = np.where((outputM <= 1.0) & (outputM > 0.9))
        
        cjInd = np.mean(Ind)
        
        if(math.isnan(cjInd)):
            cjInd = 0
        elif(cjInd < 0):
            cjInd = 0
        else:
            cjInd = cjInd

        reactionLenCJ = outputDistance[cjInd]
        #print cjInd
        #print 'Print reaction length according CJ'
        #print reactionLenCJ
        pressureCJ = outputP[cjInd]
        #print pressureCJ
        temperaturCJ = outputT[cjInd]
        #print temperaturCJ
        
    else:
        reactionLenCJ = 0;
        
    
    #############################################################
    
    if (k == 1):
        maxx = outputInd_len_ZND*5;
    else:
        maxx = outputDistance[k];

    minx = 0;

    maxT = max(outputT[:])+0.1*min(outputT[:]);
    minT = min(outputT[:])-0.1*min(outputT[:]);
    maxP = max(outputP[:])+0.1*min(outputP[:]);
    minP = min(outputP[:])-0.1*min(outputP[:]);

    
    # %Check for Eigenvalue Detonation

    if(n == a):
        print 'Error: Maximum thermicity occurs at the end of the reaction zone'
        print 'You may have an eigenvalue detonation, your final integration length may be too short,'
        print' your mixture may be too rich/lean, or something else may be wrong'
        print ' '
        print 'Mach Number (end of reaction): ' +str(outputM[a])+ ' - if close to 1, check for eigenvalue detonation'
        outputIind_time_ZND = outputTime[a]; 
        outputInd_len_ZND = outputDistance[a]; 
        outputExo_time_ZND = 0; 
        outputExo_len_ZND = 0;  
        print 'Induction Time: ' +str(outputInd_time_ZND)+'';
        print 'Exothermic Pulse Time: ' +str(outputExo_time_ZND)+'';

    elif (n == 1):
        print 'Error: Maximum thermicity occurs at the beginning of the reaction zone'
        print 'You may have an eigenvalue detonation, your final integration length may be too short,'
        print 'your mixture may be too rich/lean, or something else may be wrong'
        print ' '
        print 'Mach Number (end of reaction): ' +str(outputM[a])+ ' - if close to 1, check for eigenvalue detonation'
                                                     
        outputInd_time_ZND = outputTime[0]; 
        outputInd_len_ZND = outputDistance[0]; 
        outputExo_time_ZND = 0; 
        outputExo_len_ZND = 0;  
        print 'Induction Time: '+str(outputInd_time_ZND)+'';
        print 'Exothermic Pulse Time: '+str(outputExo_time_ZND)+'';
    else:
        print 'There is no Error'
        max_sigmadot = max(outputThermicity); # max thermicity
        half_sigmadot_flag1 = 0;
        half_sigmadot_flag2 = 0;
        # Go into a loop to find two times when sigma_dot is half its maximum
        for j in range(1,a):
            if(half_sigmadot_flag1 == 0):
                
                if(outputThermicity[j] > 0.5* max_sigmadot):
                    
                    half_sigmadot_flag1 = 1;
                    tstep1 = j;
                
            elif(half_sigmadot_flag2 == 0):
                if(outputThermicity[j] < 0.5* max_sigmadot):
                    half_sigmadot_flag2 = 1;
                    tstep2 = j;
                else:
                    tstep2 = 0;

    #Exothermic time for ZND detonation
    if(tstep2 == 0):
        print 'Error: No pulse in the thermicity'
        print '       You may have an eigenvalue detonation, your final integration length may be too short,'
        print '       your mixture may be too rich/lean, or something else may be wrong'
        outputExo_time_ZND = 0; 
        outputExo_len_ZND = 0;  
    else:
        outputExo_time_ZND = outputTime[tstep2] - outputTime[tstep1]; 
        outputExo_len_ZND = outputDistance[tstep2] - outputDistance[tstep1];
    
    '''
    print "------------------"
    print "Induction length"
    print outputInd_len_ZND
    print "------------------"
    print 'Expothermic'
    print outputExo_len_ZND
    '''

    #print "Write file"
    
    headers = "outputTime, outputP, outputRho, outputDistance, outputT, outputU, outputThermicity, outputM, outputAf, outputG, outputWt, outputSonic"
    
    np.savetxt('ZND.out', np.c_[outputTime, outputP, outputRho, outputDistance,
                                 outputT, outputU, outputThermicity, outputM,
                                 outputAf, outputG, outputWt, outputSonic],
               header=headers)     
    
    return [outputTime, outputP, outputRho, outputDistance,
            outputT, outputSpecies, outputU, maxThermicity,
            outputM, outputAf, outputG, outputWt, outputSonic,
            outputInd_time_ZND, outputInd_len_ZND, temp_grad, reactionLenCJ]

##########################################################################    

def soundSpeedM(a):
    # Matlab  SOUNDSPEED  - Speed of sound (m/s).

    '''
    #if a.isIdealGas():
    gamma = a.cp_mass()/a.cv_mass();
    wtm = a.meanMolecularWeight();
    r = GasConstant/wtm;
    b = sqrt(gamma * r * a.temperature());
    '''
    #else:
    rho0 = a.density();
    p0 = a.pressure();
    s0 = a.entropy_mass();
    rho1 = 1.001*rho0;
    a.set(Rho = rho1, S = s0);
    p1 = a.pressure();
    dpdrho_s = (p1 - p0)/(rho1 - rho0);
    b = sqrt(dpdrho_s);
    
    return b


##########################################################################
# Function to chack NAN values

def isNaN(num):
    return num != num

##########################################################################

'''
plt_num = 1
U1 = 2000;
P1 = 100000;
T1 = 300;
#q = 'H2:0.29 O2:0.1491 N2:0.5609'
rCJ = 1
#q = 'H2:2 O2:1 N2:3.76'
q = 'H2:2 O2:2 N2:7.52'
mech = 'h2air_highT.cti'

out = znd_noshk(P1, T1, q, mech, rCJ)
#out = znd_shk(U1, P1, T1, q, mech);
print out[1][16]


plt.plot(out[0], out[1])
#plt.plot(out[1][3],out[1][1])
plt.show()
'''
