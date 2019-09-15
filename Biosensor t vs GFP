import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pylab as plot
from decimal import Decimal

#Initialize some constants
vTX1 = 18.2       #reporter transcription rate constant
vTX2 = 18.2       #repressor transcription rate constant
KTX1 = 8.5        #Michaelis-Menten constant for reporter transcription 
KTX2 = 8.5        #Michaelis-Menten constant for repressor transcription
lam_m1 = 0.08     #reporter mRNA degradation rate constant
lam_m2 = 0.08     #repressor mRNA degradation rate constant
kTL1 = 0.0076       #Reporter translation rate constant, 0.0076
kTL2 = 0.0076       #Repressor translation rate constant, 0.0076
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

def model(z,t,GR,TO,TA):
    
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
    
    #ODEs  
    dmRdt = (vTX2*(GR**2)/((KTX2**2)+(GR**2))) - lam_m2*mR                                       #Repressor mRNA Conc
    dRdt = (kTL2*TLR*(mR**3)/((KTL2**3)+(mR**3)+(mF**3))) - 2*k2R*(R**2) + 2*k_2R*R2             #Repressor Conc  
    dR2dt = k2R*(R**2) - k_2R*R2 - kr*R2*O + k_r*(TO-O) - kdr1*(A**2)*R2 + k_dr1*((TA-A)/2)      #Repressor dimer Conc               
    dAdt = -2*kdr1*(A**2)*R2 + 2*k_dr1*((TA-A)/2) - 2*kdr2*(A**2)*(TO-O) + 2*k_dr2*O*((TA-A)/2)  #Analyte Conc
    dOdt = -kr*R2*O + k_r*(TO-O) + kdr2*(A**2)*(TO-O) - k_dr2*O*((TA-A)/2)                       #Operator Conc
    dmFdt = (vTX1*(O**2)/((KTX1**2)+(O**2))) - lam_m1*mF + kleak*(TO-O)                          #Reporter mRNA Conc
    dFindt = (kTL1*TLR*(mF**3)/((KTL1**3)+(mF**3)+(mR**3))) - kmat*Fin                           #Inactive Reporter Conc
    dFdt = kmat*Fin                                                                              #Active Reporter Conc  
    dTLRdt = -vlamTLR*TLR/((KlamTLR)+(TLR))                                                      #Translation Resources Conc
    
    return [dmRdt, dRdt, dR2dt, dAdt, dOdt, dmFdt, dFindt, dFdt, dTLRdt]


#Initial Conditions, conc in nM
mR0 = 0
R0 = 0
R20 = 0
R21 = 300
R2O0 = 0
R2A20 = 0
mF0 = 0
Fin0 = 0
f0 = 0
TLR0 = 1520

O0 = 8

GR0 = 8
GR3 = 0

A0 = 0                 #Lead Concentration
A1 = 1e5

TO0 = O0 + R2O0        #LOOK AT THIS----------------------<<<<<<

t = np.linspace(0,300,100000)

#-------------VARYING CONC OF LEAD----------------

F = []
F_in = []
F_O = []
F_R2 = []
F_R2O = []
F_R2A2 = []
LD = []
AZ = []
Ld = []
GFP = []
k = 0
for i in range(-4,4):
    
    j = 10**(i)
    
    TA = j + 2*R2A20
    
    z = [mR0,R0,R20,j,O0,mF0,Fin0,f0,TLR0]
    
    y,az = odeint(model,z,t,args=(GR0,TO0,TA), full_output=True, atol = 1e-21,  mxstep=5000000) #, hmax = ?? TRY Stuff....
    AZ.append(az)
    
    F.append(y[:,7])
    F_in.append(y[:,6])
    
    Ld.append(j/1000)         #Total Lead in the System
    F_R2.append(y[99999][2])  #Repressor Dimer Concentration
    LD.append(y[99999][3])    #Free Lead Concentration 
    F_O.append(y[99999][4])   #Free Operator Concentration 
    
    F_R2A2.append((TA - LD[k])/2)
    F_R2O.append(TO0 - F_O[k])
    
    plt.plot(t,F[k]*0.6022,'g-', linewidth = 0.2*(k+1))
    
    r = '%.1E' % Decimal(j)
    print('At Lead Conc ',r,', GFP : ',F[k][99999])
    
    GFP.append(F[k][99999])
    
    k+=1
    
#print(np.shape(F_R2))

    
plt.xlabel('Time (min)')
plt.ylabel('molecules of RFP') 
plt.show()

q=[]
q[:] = [A*0.6022 for A in GFP]
plt.plot(Ld,q,'bo')
plt.xscale('log')
plt.xlabel('Lead (μM)')
plt.ylabel('molecules of RFP') 
plt.show()

plt.plot(Ld,F_R2,'bo', label = 'Free Repressor Dimer')
plt.xscale('log')
plt.xlabel('Lead (μM)')
plt.ylabel('Free Repressor Dimer (nM)') 
plt.legend()
plt.show()

plt.plot(Ld,F_R2O,'go', label = 'Repressor bound to Promoter')
plt.xlabel('Lead (μM)')
#plt.xscale('log')
plt.ylabel('Repressor bound to Promoter (nM)') 
plt.legend()
plt.show()

plt.plot(Ld,F_R2A2,'ro', label = 'Repressor bound to Lead')
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Lead (μM)')
plt.ylabel('Repressor bound to Lead (nM)') 
plt.legend()
plt.show()

plt.plot(Ld,LD,'bo', label = 'Free Lead')
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Lead (μM)')
plt.ylabel('Free Lead (nM)') 
plt.legend()
plt.show()
