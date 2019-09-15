import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pylab as plot
from decimal import Decimal
from scipy.integrate import solve_ivp

def safe_division(x, y):
    if y == 0:
        return 0
    else:
        return x / y

#Initialize some constants
vTX1 = 18.2       #reporter transcription rate constant
vTX2 = 18.2       #repressor transcription rate constant
KTX1 = 8.5        #Michaelis-Menten constant for reporter transcription 
KTX2 = 8.5        #Michaelis-Menten constant for repressor transcription
lam_m1 = 0.08     #reporter mRNA degradation rate constant
lam_m2 = 0.08     #repressor mRNA degradation rate constant
vTL = 16.1       #Reporter translation rate constant, 0.0076
KTL1 = 29.9       #Michaelis-Menten constant for translation of reporter
KTL2 = 29.9       #Michaelis-Menten constant for translation of repressor
vlamTLR = 13.5    #Translation resources degradation rate constant
KlamTLR = 53.2    #Michaelis-Menten constant for degradation of TL resources
kmat = 0.244        #reporter maturation rate constant
k2R = 50          #repressor dimerization rate constant
k_2R  = 0.001     #repressor dimer dissociation rate constant
kr = 960          #association rate constant for repression
k_r = 2.4         #dissociation rate constant for repression
kdr1 = 3.0e7      #association rate constant for first derepression mechanism 
k_dr1 = 12        #dissociation rate constant for first derepression mechanism
kdr2 = 3.0e7      #association rate constant for second derepression mechanism   
k_dr2 = 4800      #dissociation rate constant for second derepression mechanism
kleak = 0.0033    #leak reporter transcription rate constant

def model(t,z):
    
    #Initialize variables in a vector
    mR = z[0]         #repressor mRNA
    R = z[1]          #repressor monomer
    R2 = z[2]         #repressor dimer               #R2O - free operator of reporter gene
    A = z[3]          #analyte                       #R2A2 - repressor-operator complex
    O = z[4]          #inactive repressor
    mF = z[5]         #reporter mRNA
    Fin = z[6]        #inactive reporter
    F = z[7]          #reporter
    TLR = z[8]        #Translation resources
    TA = z[9]
    vTL = z[10]
    
    #ODEs  
    dmRdt = (vTX2*(GR**2)/((KTX2**2)+(GR**2))) - lam_m2*mR                                       #Repressor mRNA Conc
    dRdt = (vTL*TLR*(mR**3)/((KTL2**3)+(mR**3)+(mF**3))) - 2*k2R*(R**2) + 2*k_2R*R2             #Repressor Conc  
    dR2dt = k2R*(R**2) - k_2R*R2 - kr*R2*O + k_r*(TO-O) - kdr1*(A**2)*R2 + k_dr1*((TA-A)/2)      #Repressor dimer Conc               
    dAdt = -2*kdr1*(A**2)*R2 + 2*k_dr1*((TA-A)/2) - 2*kdr2*(A**2)*(TO-O) + 2*k_dr2*O*((TA-A)/2)  #Analyte Conc
    dOdt = -kr*R2*O + k_r*(TO-O) + kdr2*(A**2)*(TO-O) - k_dr2*O*((TA-A)/2)                       #Operator Conc
    dmFdt = (vTX1*(O**2)/((KTX1**2)+(O**2))) - lam_m1*mF + kleak*(TO-O)                          #Reporter mRNA Conc
    dFindt = (vTL*TLR*(mF**3)/((KTL1**3)+(mF**3)+(mR**3))) - kmat*Fin                           #Inactive Reporter Conc
    dFdt = kmat*Fin                                                                              #Active Reporter Conc  
    dTLRdt = -vlamTLR*TLR/((KlamTLR)+(TLR))                                                      #Translation Resources Conc
    dTAdt = 0
    dvTLdt = 0
    
    return [dmRdt, dRdt, dR2dt, dAdt, dOdt, dmFdt, dFindt, dFdt, dTLRdt, dTAdt, dvTLdt]


#Initial Conditions, conc in nM
mR0 = 0
R0 = 0
R20 = 0
R21 = 300
R2O0 = 0
R2A20 = 0
mF0 = 0
Fin0 = 0
F0 = 0
TLR0 = 1520

O0 = 8
GR = 8

A0 = 0                 #Lead Concentration
A1 = 1e5

TO = O0 + R2O0        #LOOK AT THIS----------------------<<<<<<

t = np.linspace(0,300,100000)

AZ = []
LD = []

def rang(v):
    
    GFP = []
    
    k = 0
    #For loop to find Range of detection (GFP Saturation)
    for i in range(-5,30):
        
        gfp = []
        A = 10**(i/5)                          #10^(0.2) = 1.585, Expecting atleast a linear increase in fluoroscence
        TA = A + 2*R2A20
    
        z = [mR0,R0,R20,A,O0,mF0,Fin0,F0,TLR0,TA,v]
    
        #y,az = odeint(model,z,t,args=(GR0,TO0,TA,v), full_output = True, atol = 1e-20, mxstep=5000000) #, hmax = ?? TRY Stuff....
        #AZ.append(az)
        
        ans = solve_ivp(fun=model, t_span=(0,300), y0=z, method = 'BDF') #t_eval = np.linspace(0,300,100000))

        gfp = ans['y'][7]       
        GFP.append(gfp[-1])
        LD.append(A)
        
        l=0
        u=1
        
        k+=1

    
    n=0; m=0
    for n in range(1,k):
        #Check Lower Bound of Detection
        if safe_division(GFP[n],GFP[n-1]) > 1.5:
            l = GFP[n-1]
            break

    for m in range(1,k-1):
        #Check Upper Bound of Detection
        if safe_division(GFP[m],GFP[m-1]) > 1.5:
            u = GFP[m+1]
    
    return l,u



for i in range(-4,3):
    v = 0.76*10**i
    l,u = rang(v)
    r = [v, l, u]
    print(r)
    
