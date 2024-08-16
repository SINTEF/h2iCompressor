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

# Plotting files and functions:
import settingsAir as s
from plotFlow import plotFlowConditions
from plotSystem import plotSystemVariables
from plotCompressor import plotCompressorParam
from findFlow import inducerFlowWRTminimumRelativeVelocity
from pressureTest import pressureOverUnderEstimate
from plotPressure import pressurePlot4GIF

plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams.update({'font.size': 15})

Ndes = s.N0                                   # Rotational speed [rpm]     
etaStage = s.etaStage0                        # Isentropic Stage efficiency [-]           


### Inducer calculation----------------------------------------------------------------------------------
"""
For each value of Cm1 in s.Cm1i, we calculate the corresponding values of C1, T1, M1, P1, rho1, A1, rt1, U1t, W1t, Cm1(?) and beta1.
The value of Cm1 will affect the values of U1t and Ctheta1 (through ex. rt1), which in turn affects the value of W1t.
By plotting W1t against Cm1, we then find the value of Cm1 that minimises W1t.

rt1 inlet tip radius ???

"""
Ctheta1i = s.Cm1i * math.tan(math.radians(s.alpha1))            # Inducer absolute tangential velocity [degrees]          Trigonometry
C1i = (Ctheta1i ** 2 + s.Cm1i ** 2) ** 0.5                       # Inducer Absolute velocity C1 [m/s]                      Pythagoras theorem from velocity triangle 
T1i = s.T00i - (C1i ** 2) / (2 * s.Cp)                          # Inducer temperature [K]                                 Stagnation temperature relation             
M1i = C1i / ((s.k * s.R * T1i) ** 0.5)                          # Inducer Mach number [-]                                 by definition
P1i = s.P00 * (T1i / s.T00) ** (s.k / (s.k - 1))              # Inducer pressure [Pa]                                   Isentropic relation                        
rho1i = P1i / (s.R * T1i)                                        # Inducer density  [kg/m^3]                               Assuming ideal gas                          
A1i = s.mdot / (rho1i * s.Cm1i * (1 - s.B1))                   # Inducer flow area [m^2]                                 Continuity                            
#rt1i = (A1i / (math.pi * (1 - (hubtip) ** 2))) ** 0.5            # Inducer tip radius [m]
#rt1i = (A1i / math.pi + rh1 ** 2) ** 0.5                         # Inducer tip radius [m]                                  Geometry
rt1i = (A1i / (math.pi* (1 - s.rhDivr1**2) )) ** 0.5             # Inducer tip radius [m]


def inducerFlowWRTminimumRelativeVelocity(N):
    U1ti = 2 * math.pi * rt1i * N / 60                                           # Inducer blade tip speed [m/s]
    W1ti = (s.Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5                         # Inducer relative velocity [m/s]
    indice_min = np.argmin(W1ti)                                                 # Index of smallest relative velocity [m/s]
    rt1 = rt1i[indice_min]                                                       # Inducer tip radius [m]
    Ctheta1 = Ctheta1i[indice_min]                                               # Inducer angular velocity [m/s]
    C1 = C1i[indice_min]                                                         # Inducer velocity [m/s]
    T1 = T1i[indice_min]                                                         # Inducer temperature [K]
    M1 = M1i[indice_min]                                                         # Inducer Mach number [-]
    P1 = P1i[indice_min]                                                         # Inducer pressure [Pa]
    rho1 = rho1i[indice_min]                                                     # Inducer density [kg/m3]
    A1 = A1i[indice_min]                                                         # Inducer flow area [m3]
    U1t = U1ti[indice_min]                                                       # Inducer tip speed [m/s]
    Cm1 = s.Cm1i[indice_min]                                                     # Inducer meridonal velocity [m/s]
    beta1 = math.degrees(math.atan((U1t - Ctheta1) / Cm1))                       # Inducer relative velocity angle [deg]
    W1t = W1ti[indice_min]                                                       # Inducer relative velocity [m/s]
    omega = U1t/rt1                                                              # Angular velocity [rad/s]          =(2*np.pi*N/60)
    return U1ti,W1ti,rt1,Ctheta1,C1,T1,M1,P1,rho1,A1,U1t,Cm1,beta1,W1t,omega


### Impeller calculation---------------------------------------------------------------------------------
"""The isentropic enthalpy demand is used to find the approximate required work. Furthermore the  Work 
is used to find the outlet stagnation temperature which in turn is used to find the outlet stagnation pressure 
from isentropic relations. Further, geometrical relations are applied. - Petter """

# Defining iteration variables to keep track of looping


"""If the pressure estimate Prest initially over-estimates the pressure ratio by more than 5%, then the the error will only grow for increasing efficiency
        since the pressure ratio increase for increasing efficiency. Hence its necessary to check if the pressure estimate is an 
            under- or over-estimate before adjusting the efficiency eta. From then the efficiency can be either increased or decreased. """


# -----------------------  Iteration control initialization -----------------------------------


# ------------------- Initializing arrays and matrices for iteration and later avaluation -------------------
ZBarr = np.arange(s.bladeMin, s.bladeMax + 1, 1)                                  # Blade number          
beta2bArr = np.radians(np.arange( s.betamax, s.betamin, 1))      # discharge angle
beta2bArr = beta2bArr[:, np.newaxis]                         # Flipping from row to column vector

# ZB is row vector, varies along row, constant at column
# beta2b is column vector, varies along column height, constant along row distance

rhsExpLimitMat = np.exp(-8.16* np.cos( beta2bArr)/ ZBarr)                                       
trueFalseMat = np.zeros(np.shape(rhsExpLimitMat), dtype=object)
trueFalseMat = np.array([[[False, False] for _ in range(np.shape(rhsExpLimitMat)[1])] for _ in range(np.shape(rhsExpLimitMat)[0])] )
beta2ZBmat = np.zeros(np.shape(rhsExpLimitMat), dtype=object)
etaMat = np.array([[ np.nan for _ in range(np.shape(rhsExpLimitMat)[1])] for _ in range(np.shape(rhsExpLimitMat)[0])] )
sigmaWiesnerMat = 1 - (np.sqrt(np.cos(np.radians(beta2bArr))) / (ZBarr ** 0.7))
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
pressErrorMat = np.copy(etaMat)
MachExitMat = np.copy(etaMat)
WxMat = np.copy(etaMat)
Ctheta2Mat = np.copy(etaMat)
dh0SlipCorrMAt = np.copy(etaMat)
sigmaMat = np.copy(etaMat)

for iz in range(len(ZBarr)):
    for ib in range(len(beta2bArr)):
        beta2ZBmat[ib, iz] = (ZBarr[iz], beta2bArr[ib][0])

NcritArr = np.copy(etaMat)
NArr = np.copy(etaMat)

breakCheck = False


omegaCritMultr1 = np.sqrt(3*s.impellerTensileStrength/ s.impellerDensity )/(1/s.r1Divr2)
NcritMultr1 = omegaCritMultr1 *60 /(2*math.pi)
U2Crit = np.sqrt(2*s.impellerTensileStrength/s.impellerDensity)
U2 = U2Crit



""" ------------------------------------ Beginning iteration ------------------------------------ """
print('\n')

countcrit=0

dh0s = ( (s.k * s.R * s.T00) / (s.k - 1) ) * (s.Pr **((s.k - 1)/ s.k) - 1)                     # Specific isentropic compression enthalpy, constant value [kJ/kg/K]
   
""" Iterating to find a impeller tip velocity that satisfy material constraint. Finding     """
while (U2) >= U2Crit and Ndes > 0:
    # U2omega
    # U1divr2 = 2*np.pi * s.r1Divr2 * Ndes / 60     
    # secondOrderCoeff = [mu*(U2divr2**2), -np.mean(U1divr2)*Ctheta1, -Wx]
    # r2roots = np.roots(secondOrderCoeff)
    # rt2 = max(r2roots)

    U1ti, W1ti, rt1, Ctheta1, C1, T1, M1, P1, rho1, A1, U1t, Cm1, beta1, W1t, omega = inducerFlowWRTminimumRelativeVelocity(N=Ndes)
    U2 = (U1t/s.r1Divr2)

    if (U2) > s.bladeVelUpperLimit:
        break
    if (U2) >= U2Crit:
        Ndes -= 100
    rt2 = U2/omega
    Ncrit = 60*U2Crit/(2*np.pi*rt2)
    # U2 = U2divr2*rt2


""" Iterating main functionality """
for iN in range(1): # To be iterated for rpm N

    for iz in range(len(ZBarr)):
        z = ZBarr[iz]
        debug=1         # breakpoint
        
        for ib in range(len(beta2bArr)):
            expCond = rhsExpLimitMat[ib, iz]
            sigma = sigmaWiesnerMat[ib, iz]           # slip factor
            etaStage = s.etaStage0
            trueFalse2 = False
            np.rad2deg((beta2bArr[ib])[0])
            lookForBetterEta = 0
            

            D2 = 2*rt2
            rh1 = s.rhDivr1*rt1
            NDivr1 = Ndes/rt1

            
            r1tDivr2t = rt1/(rt2)
            r1r2Mat[ib, iz] = r1tDivr2t

            trueFalse1 = False
            rhs = rhsExpLimitMat[ib, iz]
            
            if r1tDivr2t < rhs:
                trueFalse1 = True
            elif r1tDivr2t > rhs:
                trueFalse1 = True
                sigma = sigma* (1  -  (   (r1tDivr2t-rhs)/(1-rhs)   )**3   )        # Correcting if ratio of radii is freater than epsilon limit
            else: 
                trueFalse1 = False
                continue

            trueFalseMat[ib, iz][0] = trueFalse1
            Vslip = (1 - sigma)*U2                              # slip velocity

            mu = sigma * s.lambda2 / (s.lambda2 - math.tan((beta2bArr[ib])[0]))                # Work input coefficient [-]    MSG: Can formula be adjusted? tan(-x=-tan(x)) 

            Ctheta2m = mu * U2                                # Absolute tangential exit velocity [m/s]         from work coefficient
            Cm2m = Ctheta2m / s.lambda2                       # Absolute meridional exit velocity [m/s]         
            C2 = (Ctheta2m ** 2 + Cm2m ** 2) ** 0.5           # Absolute exit velocity [m/s]                    from pythagoras 

            dh0SlipCorrected =  U2*(Ctheta2m-Vslip) - U1t*Ctheta1   # Alternative for finding fluid enthalpy change, for comparison
            
            while ( (trueFalseMat[ib, iz][1] == False) and (etaStage < s.etaUpperLimit and etaStage > s.etaLowerLimit) ):
                
                Wx = dh0s / etaStage                                     # Specific work [J/kg/K]
                workError = np.abs(Wx-dh0SlipCorrected)/Wx               # Comparing the two methods of finding enthalpy change


                # ------------------- Impeller outlet calculation. Stagnation properties denoted by zero. -------------------
                T02m = s.T00 + Wx * (s.k - 1) / (s.k * s.R)                                                     # Stagnation exit temperature [K]     , from dh=cp*Dt   (s.k - 1) / (s.k * s.R) = 1/cp
                P02m = s.P00 * ((etaStage * Wx * (s.k - 1) / (s.k * s.R * s.T00)) + 1) ** (s.k / (s.k - 1))     # Exit Stagnation Pressure [Pa]       , from (... todo ...)
                
                T2m = T02m - (s.k - 1) / (2 * s.k * s.R) * (C2 **2)                                             # Exit temperature [K]                            from stagnation temperature
                P2m = P02m / ((T02m / T2m) ** (s.k / (s.k - 1)))                                                # Exit pressure [Pa]                              from stagnation pressure
                rho2m = P2m / (T2m * s.R)                                                                       # Exit density [kg/m^3]                           from ideal gas law
                A2 = s.mdot / (rho2m * Cm2m)                                                                    # Exit area [m^2]                                 from continuity
                b2 = A2 / (math.pi * D2)                                                                        # Depth of impeller exit [m]                      from geometry
                M2 = U2/np.sqrt(s.k*s.R*T02m)                                                                   # Impeller exit blade mach number 
                

                # ------------------- Diffuser Calculation -------------------
                P3 = P2m + s.CpD * (P02m - P2m)               # Diffuser exit static pressure [Pa]
                C3 = C2 / s.AR                                # Diffuser exit absolute velocity [m/s]
                P03 = P3 + 0.5 * rho2m * C3 ** 2              # Diffuser exit stagnation pressure [Pa]


                # ------------------- Overall performance -------------------
                etaiterate = ( (P03 / s.P00) ** ((s.k - 1) / s.k) - 1 ) / ( (T02m / s.T00) - 1 )        # Iterative stage efficiency [-], isothermal diffuser assumed?
                Prest = ((etaiterate * U2 ** 2 * mu) / (s.Cp * T1) + 1) ** (s.k / (s.k - 1))          # Estimate of the pressure ratio, equation is validated
                PressureTestOuterLoop = (Prest-s.Pr)/s.Pr

                """Checking if teration conditions are sustained"""
                if (etaiterate > s.etaUpperLimit or etaiterate < s.etaLowerLimit) or ( etaStage > s.etaUpperLimit or etaStage < s.etaLowerLimit ):
                    trueFalse2 = False
                elif abs(PressureTestOuterLoop) > s.iterTol:
                    checkOuterLoop = False
                    trueFalse2 = False
                else:
                    trueFalse2 = True

                trueFalseMat[ib, iz][1] = trueFalse2

                if etaStage <= s.etaLowerLimit or etaStage >= s.etaUpperLimit:
                    break
                
                """ Updating values from nan's if iteration criterion is upheld"""
                if  ( trueFalse2 == True ) and U2 < s.bladeVelUpperLimit:
                    etaMat[ib, iz] = etaStage
                    U2Mat[ib, iz] = U2
                    r2Mat[ib, iz] = rt2
                    pressErrorMat[ib, iz] = abs(PressureTestOuterLoop)
                    NcritArr[ib, iz] = Ncrit
                    MachExitMat[ib, iz] = M2
                    rh1Mat[ib, iz] = rh1
                    r1Mat[ib, iz] = rt1
                    h2Mat[ib, iz] = b2
                    r2Mat[ib, iz] = rt2
                    c2Mat[ib, iz] = C2
                    PrestMat[ib, iz] = Prest
                    W1tMat[ib, iz] = W1t
                    VslipMat[ib, iz] = Vslip
                    Ctheta2Mat[ib, iz] = Ctheta2m
                    WxMat[ib, iz] = Wx
                    dh0SlipCorrMAt[ib, iz] = dh0SlipCorrected
                    sigmaMat[ib, iz] = sigma



                """ Updating efficiency for next iteration"""
                etaStage = pressureOverUnderEstimate( PressureTestOuterLoop, etaStage)
        
        
                debug = 1       # Breakpoint
            debug =1            # Breakpoint
            



# ------------- PLOTTING and DATA PRESENTATION ------------ 
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
        r"Desired $PR^*$: " + str(s.Pr) +  "\n" \
        r"$\eta _{max}$ = " + str(round(maxEta, 3)) + "   " \
        r"$\eta _{min}$ = " + str(round(minEta, 3)) + "\n" \
        r"$W_{1t}= $" + str(round(W1t, 3)) + "   " \
        r"$C_{m1}= $" + str(round(Cm1, 3)) + "\n" \
        r"$U_{1t}= $" + str(round(U1t, 3)) + "   " \
        r"$M_{1f}= $" + str(round(M1, 3)) + "\n" \
        r"$U_{2t}$ = "  +str(round(U2, 3)) + "   " \
        r"$U_{2t,crit}= $" + str(round(U2Crit, 3)) + "\n" \
         + "\n" + effText

text2 = "\n" \
        r"$\frac{r_h}{r_1}{*}$ = "  +str(s.rhDivr1) + r" $\frac{m}{m} $" +"    " \
        r"$\frac{r_1}{r_2}{*}$ = " + str(s.r1Divr2) + r" $\frac{m}{m} $" +"\n" \
        r"$r_{t1}$ = " +str(round(rt1, 3)) + r"$m$" +"   " \
        r"$r_{t2}$ = " +str(round(rt2, 3)) + r"$m$" +" \n" \
        r"$r_{h1}$ = " +str(round(rh1, 3)) + r"$m$" +" \n" \
        r"$N^{*}$ = " +str(round(Ndes, 1)) + r"$rpm$" +" \n" \
        r"$N_{crit}$ = " +str(round(Ncrit, 3)) + r"$rpm$" +" \n" \
        # r"$b_{1}$ = " +str(round(b1, 3)) + r"$m$" +" \n" \
        

""" Preparing input to plottting functions"""
systemVar = [s.etaStage0,s.lambda2, s.iterTol, ZBarr, beta2bArr]
Zflow = [Ctheta2Mat, pressErrorMat, etaMat, MachExitMat, PrestMat, trueFalseMat, WxMat, dh0SlipCorrMAt]
designParam = [rt1, rt2, rh1, Ncrit, Ndes]
flowVar = [s.Pr, W1t, Cm1, U1t, U2, U2Crit]

Zcompressor = [NcritArr, rh1Mat, h2Mat, r1Mat, r2Mat, etaMat, trueFalseMat, VslipMat]

Zsystem = [WxMat, dh0SlipCorrMAt, sigmaMat, VslipMat]
text = [text1, text2]

""" Plot of relative velocities """
ax11 = plt.figure('W1 plot')
plt.plot(s.Cm1i, W1ti, '--g')
plt.xlabel("Cm1 (m/s)")
plt.ylabel("W1t (m/s)")
plt.title("Minimisation of W1t")
plt.grid()


plotCompressorParam(systemVar, Zcompressor, designParam, flowVar, text )
plotFlowConditions(systemVar, Zflow, designParam, flowVar, text)
plotSystemVariables(systemVar, Zsystem, designParam, flowVar, text)



# plt.show(block=True)
# plt.show()

