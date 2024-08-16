"""
The following Python code takes a set of variables as inputs and calculates the geometry of an inducer, 
impeller and diffuser for a centrifugal compressor.

Reference: Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244
Author: Petter Resell (SINTEF Energy Research, 2024)

"""

### Import-------------------------------------------------------------------------------------------
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Plotting files:
from plotFlow import plotFlowConditions
from plotCompressor import plotCompressorParam

plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams.update({'font.size': 15})

### Variables-------------------------------------------------------------------------------------------
# Inlet flow parameters

mdot = 25                         # Mass flow rate [kg/s]   , mdot = ( (300* 10**3 )/( 60**2 ) )/6  

# Ndes = 120000                           # Rotational speed [rpm] 
# Ndes = 60000                          # Rotational speed [rpm]   
Ndes = 20000                          # Rotational speed [rpm]     



Cm1i = np.arange(10, 600.5, 0.5)       # Inlet Absolute meridional velocity [m/s]ø

# Fluid properties  
CpAir = 1006                            # Cp air [J/kg/K], 
MolarMassAir = 0.02897                  # Molecular weight of air [kg/mol]. 
CpH2 = 14310                            # Engineering toolbox Cp hydrogen, [J/kg/K]                                       
MolarMassH2 = 2.01568 * 10**-3          # Molecular weight of H2 [kg/mol]. 
Cp = CpH2                               # Engineering toolbox Cp hydrogen, [J/kg/K]                                       
MolarMass = MolarMassH2                 # Molecular weight of H2 [kg/mol]. 
k = 1.41                                # similar for air and hydrogen
R_uni = 8.314                           # Universal gas constant [J/mol/K]
R = R_uni / MolarMass                   # Specific gas constant [J/kg /K]

# Ti-6Al-2Sn-4Zr-2Mo titanium  properties https://www.azom.com/article.aspx?ArticleID=9298
TensileStrength = 900* 10**6       # Approximate ultimate tensile strength UTS [Pa]
impellerDensity = 1540                     # [kg/m3]

# Impeller exit parameters
lambda2 = 2                 # Exit swirl parameter                                          , check out
lambda20 = lambda2          # Exit swirl parameter, used for iteration
etaStage = 0.6             # Isentropic Stage efficiency [-]           
etaStage0 = etaStage                               

# Generall compressor parameters
# inlet conditions
P00 = 30 * (10**2) * (10**3)                    # Inlet stagnation pressure [Pa]
T00 = 293                                       # Inlet stagnation temperature [K]
Pr = 1.24                   # Pressure ratio [-]


""" ------- Vary these parameters ------- """

rhDivr1=0.55                                 # Ratio commonly given
# r1Dr2 = np.arange(0.8, 0.3, -0.05)        # rt1 divided by rt2 ,making rt2 increase through array progression
r1Dr2 = 0.35                                # rt1 divided by rt2 ,making rt2 increase through array progression

""" ------------------------------------- """


alpha1 = 0                                      # Absolute inlet velocity angle [degrees]
B1 = 0.04                                       # Boundary layer blockage [-],        substance viscocity dependant
B2 = 0.04                                       # Boundary layer blockage [-],        substance viscocity dependant NEED GREATER BLOCKAGE FOR OUTLET
AR = 2.5                                        # inlet/outlet area ratio of the diffuser [-]         
T00i = np.full(len(Cm1i), 293)                  # Inlet stagnation temperature [K]

# -------------- Diffuser Calculation --------------
etad = 0.85                                 # Estimated diffuser efficiency
CpDi = 1 - 1 / (AR ** 2)                    # Ideal pressure recovery coefficient
CpD = etad * CpDi    

### Inducer calculation----------------------------------------------------------------------------------
"""
For each value of Cm1 in Cm1i, we calculate the corresponding values of C1, T1, M1, P1, rho1, A1, rt1, U1t, W1t, Cm1(?) and beta1.
The value of Cm1 will affect the values of U1t and Ctheta1 (through ex. rt1), which in turn affects the value of W1t.
By plotting W1t against Cm1, we then find the value of Cm1 that minimises W1t.

rt1 inlet tip radius ???

"""
                                                                     # what                                       # from
Ctheta1i = Cm1i * math.tan(math.radians(alpha1))            # Inlet absolute tangential velocity [degrees]   Trigonometry
C1i = (Ctheta1i ** 2 + Cm1i ** 2) ** 0.5                    # Absolute velocity C1 [m/s]                     Pythagoras theorem from velocity triangle 
T1i = T00i - (C1i ** 2) / (2 * Cp)                          # Inlet temperature [K]                          Stagnation temperature relation             
M1i = C1i / ((k * R * T1i) ** 0.5)                          # Inlet Mach number [-]                          by definition
P1i = P00 * (T1i / T00) ** (k / (k - 1))                    # Inlet pressure [Pa]                            Isentropic relation                        
rho1i = P1i / (R * T1i)                                     # Inlet density of air [kg/m^3]                  Assuming ideal gas                          
A1i = mdot / (rho1i * Cm1i * (1 - B1))                      # Inlet flow area [m^2]                          Continuity                            
#rt1i = (A1i / (math.pi * (1 - (hubtip) ** 2))) ** 0.5      # Inlet tip radius [m]
#rt1i = (A1i / math.pi + rh1 ** 2) ** 0.5                    # Inlet tip radius [m]                           eqn (52) rapport
rt1i = (A1i / (math.pi* (1 - rhDivr1**2) )) ** 0.5 



def findRelativeInletVelocity(N):
    U1ti = 2 * math.pi * rt1i * N / 60                       # Inlet blade tip speed [m/s]                    eqn (60) rapport
    W1ti = (Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5        # Todo    : enig i denne linjen, ikke den over
    indice_min = np.argmin(W1ti)
    rt1 = rt1i[indice_min]           # Replaced by prior loop 
    Ctheta1 = Ctheta1i[indice_min]     # Get the value of Ctheta1 for minimal value of W1t
    C1 = C1i[indice_min]        
    T1 = T1i[indice_min]
    M1 = M1i[indice_min]
    P1 = P1i[indice_min]
    rho1 = rho1i[indice_min]
    A1 = A1i[indice_min]
    U1t = U1ti[indice_min]
    Cm1 = Cm1i[indice_min]
    beta1 = math.degrees(math.atan((U1t - Ctheta1) / Cm1))      # Inlet relative velocity angle [deg]
    W1t = W1ti[indice_min]
    omega = U1t/rt1                 # =2*np.pi*N/60
    return U1ti,W1ti,rt1,Ctheta1,C1,T1,M1,P1,rho1,A1,U1t,Cm1,beta1,W1t,omega


U1ti, W1ti, rt1, Ctheta1, C1, T1, M1, P1, rho1, A1, U1t, Cm1, beta1, W1t, omega = findRelativeInletVelocity(N=Ndes)

# -------------------- This bulk of code only to plot the variation of the relative velocity ------------------
def W1tIterate(Cm1i, U1ti, Ctheta1i):                              
    #return (Cm1i ** 2 + ((U1ti ** 2) - (Ctheta1i ** 2))) ** 0.5     # Todo MSG: Should be (Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5?
    return (Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5


W1tIter = W1tIterate(Cm1i, U1ti, Ctheta1i)      # Get W1t as a function of Cm1

plt.figure('W1 plot')
plt.plot(Cm1i, W1tIter)
plt.xlabel("Cm1 (m/s)")
plt.ylabel("W1t (m/s)")
plt.grid()
plt.title("Minimisation of W1t")

# -------------------------------------------------------------------------------------------------------------

# Iterate-function decides which  variable to change when iterating. Change this by changing default arguments a1, a2 & a3 


### Impeller calculation---------------------------------------------------------------------------------
"""The isentropic enthalpy demand is used to find the approximate required work. Furthermore the  Work 
is used to find the outlet stagnation temperature which in turn is used to find the outlet stagnation pressure 
from isentropic relations. Further, geometrical relations are applied. - Petter """

# Defining iteration variables to keep track of looping


"""If the pressure estimate Prest initially over-estimates the pressure ratio by more than 5%, then the the error will only grow for increasing efficiency
        since the pressure ratio increase for increasing efficiency. Hence its necessary to check if the pressure estimate is an 
            under- or over-estimate before adjusting the efficiency eta. From then the efficiency can be either increased or decreased. """

def pressureOverUnderEstimate( pressureTest, etaStage):
    if pressureTest == np.nan:
        debug = 0                   # breakpoint for debugging
    if  pressureTest > 0:
        etaST = etaStage - 0.003                  # todo: need better way of iterating etaStage
    elif  pressureTest < 0:
        etaST = etaStage + 0.003

    return etaST

### -----------------------  Iteration control initialization -----------------------------------
etaLowerLimit = 0.25                        # Lowest efficiency allowed
etaUpperLimit = 0.85                        # Highest efficiency allowed
bladeVelUpperLimit = 1200                   # Highest blade velocity allowed
bladeVelLowerLimit = 0                      # lowest blade velocity allowed 
betamax = -45                               # "Maximum" beta iterated over
betamin = 0                                 # "Minimum" beta iterated over  
bladeMin = 1                                # Lowest bladenumber allowed
bladeMax = 30                               # Highest bladenumber allowed
iterTol = 0.005                            # loop tolerance condition

# ------------------- Initializing arrays and matrices for iteration and later avaluation -------------------
lambda2Arr = np.arange(0.2)                                  # Exit swirl number
ZBarr = np.arange(1, 31, 1)                                  # Blade number          
beta2bArr = np.radians(np.arange( betamax, betamin, 1))      # discharge angle
beta2bArr = beta2bArr[:, np.newaxis]                         # Flipping from row to column vector

# ZB is row vector, varies along row, constant at column
# beta2b is column vector, varies along column height, constant along row distance

rhsExpLimitMat = np.exp(-8.16* np.cos( beta2bArr)/ ZBarr)                                       
trueFalseMat = np.zeros(np.shape(rhsExpLimitMat), dtype=object)
trueFalseMat = np.array([[[False, False] for _ in range(np.shape(rhsExpLimitMat)[1])] for _ in range(np.shape(rhsExpLimitMat)[0])] )
beta2ZBmat = np.zeros(np.shape(rhsExpLimitMat), dtype=object)
etaMat = np.array([[ np.nan for _ in range(np.shape(rhsExpLimitMat)[1])] for _ in range(np.shape(rhsExpLimitMat)[0])] )
sigmaWmat = 1 - (np.sqrt(np.cos(np.radians(beta2bArr))) / (ZBarr ** 0.7))
r1r2Mat = np.zeros(np.shape(rhsExpLimitMat), dtype=object)
r2Mat = np.copy(etaMat)
r1Mat = np.copy(etaMat)
rh1Mat = np.copy(etaMat)
h2Mat = np.copy(etaMat)
U2Mat = np.copy(etaMat)
c2Mat = np.copy(etaMat)
U1Mat = np.copy(etaMat)
PrestMat = np.copy(etaMat)
W1tMat = np.copy(etaMat)
VslipMat = np.copy(etaMat)

pressErrorMat = np.copy(etaMat)

for iz in range(len(ZBarr)):
    for ib in range(len(beta2bArr)):
        beta2ZBmat[ib, iz] = (ZBarr[iz], beta2bArr[ib][0])

NcritArr = np.copy(etaMat)
NArr = np.copy(etaMat)
MachExitArr = np.copy(etaMat)

breakCheck = False


omegaCritMultr1 = np.sqrt(3*TensileStrength/ impellerDensity )/(1/r1Dr2)
NcritMultr1 = omegaCritMultr1 *60 /(2*math.pi)
U2Crit = np.sqrt(2*TensileStrength/impellerDensity)
U2 = U2Crit



""" ------------------------------------ Beginning iteration ------------------------------------ """
print('\n')

countcrit=0

dh0s = ( (k * R * T00) / (k - 1) ) * (Pr **((k - 1)/ k) - 1)                     # Specific isentropic compression enthalpy, constant value [kJ/kg/K]
   
while (U2) >= U2Crit and Ndes > 0:
    # U2omega
    # U1divr2 = 2*np.pi * r1Dr2 * Ndes / 60     
    # secondOrderCoeff = [mu*(U2divr2**2), -np.mean(U1divr2)*Ctheta1, -Wx]
    # r2roots = np.roots(secondOrderCoeff)
    # rt2 = max(r2roots)

    U1ti, W1ti, rt1, Ctheta1, C1, T1, M1, P1, rho1, A1, U1t, Cm1, beta1, W1t, omega = findRelativeInletVelocity(N=Ndes)
    U2 = (U1t/r1Dr2)

    rt2 = U2/omega
    Ntest = 60*U2/(2*np.pi*rt2)
    Ncrit = 60*U2Crit/(2*np.pi*rt2)
    # U2 = U2divr2*rt2
    if (U2) > bladeVelUpperLimit:
        break
    if (U2) >= U2Crit:
        Ndes -= 100

for iN in range(1): # To be iterated for rpm N

    for iz in range(len(ZBarr)):
        z = ZBarr[iz]
        debug=1         # breakpoint
        
        for ib in range(len(beta2bArr)):
            expCond = rhsExpLimitMat[ib, iz]
            sigma = sigmaWmat[ib, iz]           # slip factor
            etaStage = etaStage0
            trueFalse1 = False
            trueFalse2 = False
            b = np.rad2deg((beta2bArr[ib])[0])
            lookForBetterEta = 0
            

            D2 = 2*rt2
            rh1 = rhDivr1*rt1
            NDivr1 = Ndes/rt1

            
            r1tDivr2t = rt1/(rt2)
            r1r2Mat[ib, iz] = r1tDivr2t

            trueFalse1 = False
            rhs = rhsExpLimitMat[ib, iz]
            
            if r1tDivr2t < rhs:
                trueFalse1 = True
            elif r1tDivr2t > rhs:
                trueFalse1 = True
                sigma = sigma* (1  -  (   (r1tDivr2t-rhs)/(1-rhs)   )**3   )
            else: 
                trueFalse1 = False
                continue

            trueFalseMat[ib, iz][0] = trueFalse1
            Vslip = (1 - sigma)*U2

            mu = sigma * lambda2 / (lambda2 - math.tan((beta2bArr[ib])[0]))                # Work input coefficient [-]    MSG: Can formula be adjusted? tan(-x=-tan(x)) 

            Ctheta2m = mu * U2                                # Absolute tangential exit velocity [m/s]         from work coefficient
            Cm2m = Ctheta2m / lambda2                         # Absolute meridional exit velocity [m/s]         
            C2 = (Ctheta2m ** 2 + Cm2m ** 2) ** 0.5           # Absolute exit velocity [m/s]                    from pythagoras 
            
            while ( (trueFalseMat[ib, iz][1] == False) and (etaStage < etaUpperLimit and etaStage > etaLowerLimit) ):

                Wx = dh0s / etaStage                                                             # Specific work [J/kg/K]
                U222 = ((U1t * Ctheta1 + Wx) / mu) ** 0.5    # Exit Blade Speed [m/s]
                U222 = ((U1t * Ctheta1 + Wx) / Ctheta2m)    # Exit Blade Speed [m/s]
                T02m = T00 + Wx * (k - 1) / (k * R)                                              # Stagnation exit temperature [K]     , from dh=cp*Dt   (k - 1) / (k * R) = 1/cp
                P02m = P00 * ((etaStage * Wx * (k - 1) / (k * R * T00)) + 1) ** (k / (k - 1))    # Exit Stagnation Pressure [Pa]       , from (... todo ...)
                
                T2m = T02m - (k - 1) / (2 * k * R) * (C2 **2)     # Exit temperature [K]                            from stagnation temperature
                P2m = P02m / ((T02m / T2m) ** (k / (k - 1)))      # Exit pressure [Pa]                              from stagnation pressure
                rho2m = P2m / (T2m * R)                           # Exit density [kg/m^3]                           from ideal gas law
                A2 = mdot / (rho2m * Cm2m)                        # Exit area [m^2]                                 from continuity
                b2 = A2 / (math.pi * D2)                          # Depth of impeller exit [m]                      from geometry
                M2 = U2/np.sqrt(k*R*T02m)
                

                # Diffuser Calculation----------------------------------------------------------------------------------
                P3 = P2m + CpD * (P02m - P2m)               # Diffuser exit static pressure [Pa]
                C3 = C2 / AR                                # Diffuser exit absolute velocity [m/s]
                P03 = P3 + 0.5 * rho2m * C3 ** 2            # Diffuser exit stagnation pressure [Pa]


                ### Overall performance----------------------------------------------------------------------------------
                etaiterate = ( (P03 / P00) ** ((k - 1) / k) - 1 ) / ( (T02m / T00) - 1 )        # Iterative stage efficiency [-], isothermal diffuser assumed?
                Prest = ((etaiterate * U2 ** 2 * mu) / (Cp * T1) + 1) ** (k / (k - 1))          # Estimate of the pressure ratio, equation is validated
                PressureTestOuterLoop = (Prest-Pr)/Pr
                if ZBarr[iz]>20:
                    debug=0
                if (etaiterate > etaUpperLimit or etaiterate < etaLowerLimit) or ( etaStage > etaUpperLimit or etaStage < etaLowerLimit ):
                    trueFalse2 = False
                elif abs(PressureTestOuterLoop) > iterTol:
                    checkOuterLoop = False
                    trueFalse2 = False
                else:
                    trueFalse2 = True

                trueFalseMat[ib, iz][1] = trueFalse2

                if etaStage <= etaLowerLimit or etaStage >= etaUpperLimit:
                    break
                
                if ( trueFalse1 == True ) and ( trueFalse2 == True ) and U2 < bladeVelUpperLimit:
                    etaMat[ib, iz] = etaStage
                    U2Mat[ib, iz] = U2
                    r2Mat[ib, iz] = rt2
                    pressErrorMat[ib, iz] = abs(PressureTestOuterLoop)
                    NcritArr[ib, iz] = Ncrit
                    MachExitArr[ib, iz] = M2
                    rh1Mat[ib, iz] = rh1
                    r1Mat[ib, iz] = rt1
                    h2Mat[ib, iz] = b2
                    r2Mat[ib, iz] = rt2
                    c2Mat[ib, iz] = C2
                    PrestMat[ib, iz] = Prest
                    W1tMat[ib, iz] = W1t
                    VslipMat[ib, iz] = Vslip


                etaStage = pressureOverUnderEstimate( PressureTestOuterLoop, etaStage)
        
        
                debug = 1       #breakpoint
            debug =1
            



# ------------- PLOTTING ------------ 

x = (ZBarr)                     # SAME FOR ALL COUNTOURS
y = np.rad2deg((beta2bArr))     # SAME FOR ALL COUNTOURS
X, Y = np.meshgrid(x, y)  
lvls = 20
colorTheme = 'Reds'

# ---- Old plotting inside here ----
"""
fig, axs21 = plt.subplots(3, 3)
fig.set_figwidth(15)
fig.set_figheight(15)
fig.tight_layout(pad=7.0)
fig.suptitle(r'Flow properties for proposed  $\eta$ = ' + str(etaStage0) +  r' ,  $\lambda _{2}$ = ' + str(lambda2) + ' and Tolerance set to ' +str(iterTol*100) + '% .', y=0.98, x=0.38)
fig.subplots_adjust(top=0.9, bottom=0.09)



# ---------------- Flow outlet Velocity plot ---------------
i=1
j=1
Z = c2Mat
con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
cbar = fig.colorbar(con, ax=axs21[i, j])
cbar.ax.tick_params(labelsize=10)
axs21[i, j].invert_yaxis()
axs21[i, j].set_xticks(x[1::4])
axs21[i, j].set_xticklabels(x[1::4], fontsize=10)
axs21[i, j].set_yticks(np.arange(betamax, betamin , 4))
axs21[i, j].set_yticklabels(np.arange(betamax, betamin , 4), fontsize=10)
axs21[i, j].set_xlabel(r'Blade number $Z_B$ ' , fontsize=12)
axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
axs21[i, j].set_title('Impeller outlet flow velocity [m/s]', fontsize=12)
axs21[i, j].grid()


# ---------------- Pressure estimate error plot ---------------
i=0 
j=2
Z = pressErrorMat
con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
cbar = fig.colorbar(con, ax=axs21[i, j])
cbar.ax.tick_params(labelsize=10)
axs21[i, j].invert_yaxis()
axs21[i, j].set_xticks(x[1::4])
axs21[i, j].set_xticklabels(x[1::4], fontsize=10)
axs21[i, j].set_yticks(np.arange(betamax, betamin, 4))
axs21[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs21[i, j].set_xlabel(r'Blade number $Z_B$ [deg]' , fontsize=12)
axs21[i, j].set_ylabel(r'$ \beta _{2B}$', fontsize=12)
axs21[i, j].set_title(r'Pressure estimate error contour plot [-]', fontsize=12)
axs21[i, j].grid()

# -------------- Normalized Efficiency plot -------------
i=0
j=1
Z = etaMat/etaStage                  
con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
cbar = fig.colorbar(con, ax=axs21[i, j])
cbar.ax.tick_params(labelsize=10)
axs21[i, j].invert_yaxis()
axs21[i, j].set_xticks(x[1::4])
axs21[i, j].set_xticklabels(x[1::4], fontsize=10)
axs21[i, j].set_yticks(np.arange(betamax, betamin, 4))
axs21[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
axs21[i, j].set_title('Normalized required efficiency [-]' , fontsize=12)
axs21[i, j].grid()

# -------------- Efficiency plot -------------
i=0
j=0
Z = etaMat                  
con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
cbar = fig.colorbar(con, ax=axs21[i, j])
cbar.ax.tick_params(labelsize=10)
axs21[i, j].invert_yaxis()
axs21[i, j].set_xticks(x[1::4])
axs21[i, j].set_xticklabels(x[1::4], fontsize=10)
axs21[i, j].set_yticks(np.arange(betamax, betamin, 4)) #, fontsize=10)
axs21[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
axs21[i, j].set_ylabel(r' $ \beta _{2B}$ [deg]', fontsize=12)
axs21[i, j].set_title(r'Required efficiency [-]', fontsize=12)
axs21[i, j].grid()

# -------------- Mach plot -------------
i=1
j=0
Z = MachExitArr                  
con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
cbar = fig.colorbar(con, ax=axs21[i, j])
cbar.ax.tick_params(labelsize=10)
axs21[i, j].invert_yaxis()
axs21[i, j].set_xticks(x[1::4])
axs21[i, j].set_xticklabels(x[1::4], fontsize=10)
axs21[i, j].set_yticks(np.arange(betamax, betamin, 4)) #, fontsize=10)
axs21[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
axs21[i, j].set_title('Exit blade Mach number [-] ' , fontsize=12)
axs21[i, j].grid()

# -------------- PR plot -------------
i=1
j=2
Z = PrestMat                  
con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
cbar = fig.colorbar(con, ax=axs21[i, j])
cbar.ax.tick_params(labelsize=10)
axs21[i, j].invert_yaxis()
axs21[i, j].set_xticks(x[1::4])
axs21[i, j].set_xticklabels(x[1::4], fontsize=10)
axs21[i, j].set_yticks(np.arange(betamax, betamin, 4)) #, fontsize=10)
axs21[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
axs21[i, j].set_title('Estimated PR [-] ' , fontsize=12)
axs21[i, j].grid()
# -------------- PR plot -------------
i=2
j=0
Z = VslipMat                  
con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
cbar = fig.colorbar(con, ax=axs21[i, j])
cbar.ax.tick_params(labelsize=10)
axs21[i, j].invert_yaxis()
axs21[i, j].set_xticks(x[1::4])
axs21[i, j].set_xticklabels(x[1::4], fontsize=10)
axs21[i, j].set_yticks(np.arange(betamax, betamin, 4)) #, fontsize=10)
axs21[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
axs21[i, j].set_title('Slip velocity  ' , fontsize=12)
axs21[i, j].grid()
# -------------- Textbox -------------
i=2
j=1
maxEta = np.nanmax(etaMat)
minEta = np.nanmin(etaMat)
countTrue = np.sum(np.all(trueFalseMat, axis=-1))
totalCases = len(ZBarr)*len(beta2bArr)
text1 = r"Valid cases: " + str(countTrue) + '/'  +str(totalCases) + "\n" \
        r"Desired PR: " + str(Pr) +  "\n" \
        r"$\eta _{max}$ = " + str(round(maxEta, 3)) + "   " \
        r"$\eta _{min}$ = " + str(round(minEta, 3)) + "\n" \
        r"$W_{1t}= $" + str(round(W1t, 3)) + "   " \
        r"$C_{m1}= $" + str(round(Cm1, 3)) + "\n" \
        r"$U_{1t}= $" + str(round(U1t, 3)) + "\n" \
        r"$U2$ = "  +str(round(U2, 3)) + "   " \
        r"$U_{2,crit}= $" + str(round(U2Crit, 3)) + "\n"
axs21[i, j].axis('off')  # Turn off the axis
axs21[i, j].text(0.0, 0.5, text1, ha='left', va='center', fontsize=12, linespacing = 1.8 )

# -------------- Textbox -------------
i=2
j=2
text2 = r"$\frac{r_h}{r_1}$ = "  +str(rhDivr1) + r" $\frac{m}{m} $" +"    " \
        r"$\frac{r_1}{r_2}$ = " + str(r1Dr2) + r" $\frac{m}{m} $" +"\n" \
        r"$r_{t1}$ = " +str(round(rt1, 3)) + r"$m$" +" \n" \
        r"$r_{h1}$ = " +str(round(rh1, 3)) + r"$m$" +" \n" \
        r"$r_{t2}$ = " +str(round(rt2, 3)) + r"$m$" +" \n" \
        r"$N_{crit}$ = " +str(round(Ncrit, 3)) + r"$m$" +" \n" \
        
axs21[i, j].axis('off')  # Turn off the axis
axs21[i, j].text(0.0, 0.5, text2, ha='left', va='center', fontsize=12, linespacing =1.8 )


#----------------------------------------------------------------------
#---------------------------- New Figure ------------------------------
#----------------------------------------------------------------------


fig, axs2 = plt.subplots(3, 3)
fig.set_figwidth(15)
fig.set_figheight(15)
fig.suptitle(r'Compressor design parameters for proposed $\eta$ = ' + str(etaStage0) +  r' , $\lambda _{2}$ = ' + str(lambda2) + ' and Tolerance set to ' +str(iterTol*100) + '% .', y=0.98, x=0.4)
fig.tight_layout(pad=7.0)
fig.subplots_adjust(top=0.9, bottom=0.09)


# -------------- RPM plot -------------
i=0
j=2
Z = NcritArr                  
con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
cbar = fig.colorbar(con, ax=axs2[i, j])
cbar.ax.tick_params(labelsize=10)
axs2[i, j].invert_yaxis()
axs2[i, j].set_xticks(x[1::4])
axs2[i, j].set_xticklabels(x[1::4], fontsize=10)
axs2[i, j].set_yticks(np.arange(betamax, betamin, 4)) #, fontsize=10)
axs2[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
axs2[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
axs2[i, j].set_title(r'Critical RPM ' , fontsize=12)
axs2[i, j].grid()



# -------------- Inlet hub radius plot -------------
i=0
j=1
Z = rh1Mat                  
con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
cbar = fig.colorbar(con, ax=axs2[i, j])
cbar.ax.tick_params(labelsize=10)
axs2[i, j].invert_yaxis()
axs2[i, j].set_xticks(x[1::4])
axs2[i, j].set_xticklabels(x[1::4], fontsize=10)
axs2[i, j].set_yticks(np.arange(betamax, betamin, 4)) #, fontsize=10)
axs2[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
axs2[i, j].set_ylabel(r' $ \beta _{2B}$ [deg]', fontsize=12)
axs2[i, j].set_title('Inlet hub radius [m]' , fontsize=12)
axs2[i, j].grid()

# -------------- Outlet exit height plot -------------
i=1
j=1
Z = h2Mat                  
con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
cbar = fig.colorbar(con, ax=axs2[i, j])
cbar.ax.tick_params(labelsize=10)
axs2[i, j].invert_yaxis()
axs2[i, j].set_xticks(x[1::4])
axs2[i, j].set_xticklabels(x[1::4], fontsize=10)
axs2[i, j].set_yticks(np.arange(betamax, betamin, 4)) #, fontsize=10)
axs2[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
axs2[i, j].set_ylabel(r' $ \beta _{2B}$ [deg]', fontsize=12)
axs2[i, j].set_title('Outlet area height [m]', fontsize=12)
axs2[i, j].grid()

# -------------- inducer radius plot -------------
i=0
j=0
Z = r1Mat              
con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
cbar = fig.colorbar(con, ax=axs2[i, j])
cbar.ax.tick_params(labelsize=10)
axs2[i, j].invert_yaxis()
axs2[i, j].set_xticks(x[1::4])
axs2[i, j].set_xticklabels(x[1::4], fontsize=10)
axs2[i, j].set_yticks(np.arange(betamax, betamin, 4)) #, fontsize=10)
axs2[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
axs2[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
axs2[i, j].set_title('Inducer radius [m]' , fontsize=12)
axs2[i, j].grid()

# ---------------- Impeller exit Radius plot ---------------
i=1
j=0
Z = r2Mat
con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
cbar = fig.colorbar(con, ax=axs2[i, j])
cbar.ax.tick_params(labelsize=10)
axs2[i, j].invert_yaxis()
axs2[i, j].set_xticks(x[1::4])
axs2[i, j].set_xticklabels(x[1::4], fontsize=10)
axs2[i, j].set_yticks(np.arange(betamax, betamin, 4))
axs2[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs2[i, j].set_xlabel(r'Blade number $Z_B$ [deg]', fontsize=12 )
axs2[i, j].set_ylabel(r' $ \beta _{2B}$', fontsize=12)
axs2[i, j].set_title(r'Impeller exit radius [m]', fontsize=12)
axs2[i, j].grid()

# -------------- Textbox -------------
i=2
j=1
maxEta = np.nanmax(etaMat)
minEta = np.nanmin(etaMat)
countTrue = np.sum(np.all(trueFalseMat, axis=-1))
totalCases = len(ZBarr)*len(beta2bArr)
text1 = r"Valid cases: " + str(countTrue) + '/'  +str(totalCases) + "\n" \
        r"Desired PR: " + str(Pr) +  "\n" \
        r"$\eta _{max}$ = " + str(round(maxEta, 3)) + "   " \
        r"$\eta _{min}$ = " + str(round(minEta, 3)) + "\n" \
        r"$W_{1t}= $" + str(round(W1t, 3)) + "   " \
        r"$C_{m1}= $" + str(round(Cm1, 3)) + "\n" \
        r"$U_{1t}= $" + str(round(U1t, 3)) + "\n" \
        r"$U2$ = "  +str(round(U2, 3)) + "   " \
        r"$U_{2,crit}= $" + str(round(U2Crit, 3)) + "\n"

text2 = r"$\frac{r_h}{r_1}$ = "  +str(rhDivr1) + r" $\frac{m}{m} $" +"    " \
        r"$\frac{r_1}{r_2}$ = " + str(r1Dr2) + r" $\frac{m}{m} $" +"\n" \
        r"$r_{t1}$ = " +str(round(rt1, 3)) + r"$m$" +" \n" \
        r"$r_{h1}$ = " +str(round(rh1, 3)) + r"$m$" +" \n" \
        r"$r_{t2}$ = " +str(round(rt2, 3)) + r"$m$" +" \n" \
        r"$N_{crit}$ = " +str(round(Ncrit, 3)) + r"$m$" +" \n" \
        
axs2[i, j].axis('off')  # Turn off the axis
axs2[i, j].text(0.0, 0.5, text1, ha='left', va='center', fontsize=12, linespacing =1.8 )

# -------------- Textbox -------------
i=2
j=2

axs2[i, j].axis('off')  # Turn off the axis
axs2[i, j].text(0.0, 0.5, text2, ha='left', va='center', fontsize=12, linespacing =1.8 )


"""
# ---- Old plotting inside here ----


""" Can check how many cases is within a certain error margin 
for iz in range(np.shape(pressErrorMat)[1]):
    for ib in range(np.shape(pressErrorMat)[0]):
        if pressErrorMat[ib, iz] < 0.0001: 
            bestBladesAndAngles.append( (ZBarr[iz], np.rad2deg(beta2bArr[ib])) )

print(len(bestBladesAndAngles)) """


# ------ New plotting using funcitons ------
# Finding best efficiency point
minimumError = (np.nanmin((pressErrorMat)))
minimumIndice = ( np.where(pressErrorMat == minimumError) )
IBmin = int(minimumIndice[0][0])
IZmin = int(minimumIndice[1][0]) 
betaMinErr = round(np.rad2deg(float(beta2bArr[IBmin][0])), 3)
zMinErr = round(float(ZBarr[IZmin]),3)
etaMinErr = round(float(etaMat[IBmin, IZmin]), 3)
effText = r"Minimum error using $Z_{B}$ = " + str(zMinErr) + r" and $\beta _{2b}$ = " + str(betaMinErr) + r" with efficiency $\eta$ = " +str(etaMinErr) 

# making nice text to go with plots
maxEta = np.nanmax(etaMat)
minEta = np.nanmin(etaMat)
countTrue = np.sum(np.all(trueFalseMat, axis=-1))
totalCases = len(ZBarr)*len(beta2bArr)
text1 = "\n" \
        r"Valid cases: " + str(countTrue) + '/'  +str(totalCases) + "\n" \
        r"Desired PR: " + str(Pr) +  "\n" \
        r"$\eta _{max}$ = " + str(round(maxEta, 3)) + "   " \
        r"$\eta _{min}$ = " + str(round(minEta, 3)) + "\n" \
        r"$W_{1t}= $" + str(round(W1t, 3)) + "   " \
        r"$C_{m1}= $" + str(round(Cm1, 3)) + "\n" \
        r"$U_{1t}= $" + str(round(U1t, 3)) + "\n" \
        r"$U2$ = "  +str(round(U2, 3)) + "   " \
        r"$U_{2,crit}= $" + str(round(U2Crit, 3)) + "\n" \
         + "\n" + effText

text2 = "\n" \
        r"$\frac{r_h}{r_1}$ = "  +str(rhDivr1) + r" $\frac{m}{m} $" +"    " \
        r"$\frac{r_1}{r_2}$ = " + str(r1Dr2) + r" $\frac{m}{m} $" +"\n" \
        r"$r_{t1}$ = " +str(round(rt1, 3)) + r"$m$" +" \n" \
        r"$r_{h1}$ = " +str(round(rh1, 3)) + r"$m$" +" \n" \
        r"$r_{t2}$ = " +str(round(rt2, 3)) + r"$m$" +" \n" \
        r"$N$ = " +str(round(Ndes, 1)) + r"$m$" +" \n" \
        r"$N_{crit}$ = " +str(round(Ncrit, 3)) + r"$m$" +" \n" \
        
systemVar = [etaStage0, lambda2, iterTol, ZBarr, beta2bArr]
Zflow = [c2Mat, pressErrorMat, etaMat, MachExitArr, PrestMat, VslipMat, trueFalseMat]
designParam = [rt1, rt2, rh1, Ncrit, Ndes]
flowVar = [Pr, W1t, Cm1, U1t, U2, U2Crit]

Zcompressor = [NcritArr, rh1Mat, h2Mat, r1Mat, r2Mat, etaMat, trueFalseMat]

text = [text1, text2]

plotCompressorParam(systemVar, Zcompressor, designParam, flowVar, text )
plotFlowConditions(systemVar, Zflow, designParam, flowVar, text)

bestBladesAndAngles = []



plt.show(block=True)

