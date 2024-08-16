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
import settings
# import logic

# Plotting files and functions:
from plotFlow import plotFlowConditions
from plotSystem import plotSystemVariables
from plotCompressor import plotCompressorParam
from plotVelocities import plotVelocities
from plotText import plotText
# from findFlow import inducerFlowWRTminimumRelativeVelocity
from pressureTest import pressureOverUnderEstimate
from plotPressure import pressurePlot4GIF

plt.rcParams.update(plt.rcParamsDefault)            # For plotting
plt.rcParams.update({'font.size': 15})

# Ndes = settings.N0                                   # Rotational speed [rpm]     
Ndes = settings.N0                                   # Rotational speed [rpm]     
etaStage = settings.etaStage0                        # Isentropic Stage efficiency [-]           


### Inducer calculation----------------------------------------------------------------------------------
"""
For each value of Cm1 in settings.Cm1i, we calculate the corresponding values of C1, T1, M1, P1, rho1, A1, r1, U1t, W1t, Cm1(?) and beta1.
The value of Cm1 will affect the values of U1t and Ctheta1 (through ex. r1), which in turn affects the value of W1t.
By plotting W1t against Cm1, we then find the value of Cm1 that minimises W1t.

r1 inlet tip radius ???

"""
Ctheta1i = settings.Cm1i * math.tan(math.radians(settings.alpha1))            # Inducer absolute tangential velocity [degrees]          Trigonometry
C1i = (Ctheta1i ** 2 + settings.Cm1i ** 2) ** 0.5                       # Inducer Absolute velocity C1 [m/s]                      Pythagoras theorem from velocity triangle 
T1i = settings.T00i - (C1i ** 2) / (2 * settings.Cp)                          # Inducer temperature [K]                                 Stagnation temperature relation             
M1i = C1i / ((settings.k * settings.R * T1i) ** 0.5)                          # Inducer Mach number [-]                                 by definition
P1i = settings.P00 * (T1i / settings.T00) ** (settings.k / (settings.k - 1))              # Inducer pressure [Pa]                                   Isentropic relation                        
rho1i = P1i / (settings.R * T1i)                                        # Inducer density  [kg/m^3]                               Assuming ideal gas                          
A1i = settings.mdot / (rho1i * settings.Cm1i * (1 - settings.B1))                   # Inducer flow area [m^2]                                 Continuity                                                        
rt1i = (A1i / (math.pi* (1 - settings.rhDivr1**2) )) ** 0.5             # Inducer tip radius [m]



def inducerFlowWRTminimumRelativeVelocity(N):
    U1ti = 2 * math.pi * rt1i * N / 60                                           # Inducer blade tip speed [m/s]
    W1ti = (settings.Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5                         # Inducer relative velocity [m/s]
    indice_min = np.argmin(W1ti)                                                 # Index of smallest relative velocity [m/s]
    r1 = rt1i[indice_min]                                                       # Inducer tip radius [m]
    Ctheta1 = Ctheta1i[indice_min]                                               # Inducer angular velocity [m/s]
    C1 = C1i[indice_min]                                                         # Inducer velocity [m/s]
    T1 = T1i[indice_min]                                                         # Inducer temperature [K]
    M1 = M1i[indice_min]                                                         # Inducer Mach number [-]
    P1 = P1i[indice_min]                                                         # Inducer pressure [Pa]
    rho1 = rho1i[indice_min]                                                     # Inducer density [kg/m3]
    A1 = A1i[indice_min]                                                         # Inducer flow area [m3]
    U1t = U1ti[indice_min]                                                       # Inducer tip speed [m/s]
    Cm1 = settings.Cm1i[indice_min]                                                     # Inducer meridonal velocity [m/s]
    beta1 = math.degrees(math.atan((U1t - Ctheta1) / Cm1))                       # Inducer relative velocity angle [deg]
    W1t = W1ti[indice_min]                                                       # Inducer relative velocity [m/s]
    omega = U1t/r1                                                              # Angular velocity [rad/s]          =(2*np.pi*N/60)
    return U1ti,W1ti,r1,Ctheta1,C1,T1,M1,P1,rho1,A1,U1t,Cm1,beta1,W1t,omega


### Impeller calculation---------------------------------------------------------------------------------
"""The isentropic enthalpy demand is used to find the approximate required work. Furthermore the  Work 
is used to find the outlet stagnation temperature which in turn is used to find the outlet stagnation pressure 
from isentropic relationsettings. Further, geometrical relations are applied. - Petter """


"""If the pressure estimate Prest initially over-estimates the pressure ratio by more than 5%, then the the error will only grow for increasing efficiency
        since the pressure ratio increase for increasing efficiency. Hence its necessary to check if the pressure estimate is an 
            under- or over-estimate before adjusting the efficiency eta. From then the efficiency can be either increased or decreased. """


# -----------------------  Iteration control initialization -----------------------------------


# ------------------- Initializing arrays and matrices for iteration and later avaluation -------------------
ZBarr = np.arange(settings.bladeMin, settings.bladeMax + 1, 1)                                  # Blade number          
beta2BArr = np.radians(np.arange( settings.beta2Bmax, settings.beta2Bmin, 1))      # discharge angle
beta2BArr = beta2BArr[:, np.newaxis]                         # Flipping from row to column vector

# ZB is row vector, varies along row, constant at column
# beta2B is column vector, varies along column height, constant along row distance

""" Making matrices for iteration later. Filling with nans that are only replaced if all conditions are met. """
rhsExpLimitMat = np.exp(-8.16* np.cos( beta2BArr)/ ZBarr)                                       
trueFalseMat = np.zeros(np.shape(rhsExpLimitMat), dtype=object)
trueFalseMat = np.array([[[False, False] for _ in range(np.shape(rhsExpLimitMat)[1])] for _ in range(np.shape(rhsExpLimitMat)[0])] )
etaMat = np.array([[ np.nan for _ in range(np.shape(rhsExpLimitMat)[1])] for _ in range(np.shape(rhsExpLimitMat)[0])] )
sigmaWiesnerMat = 1 - (np.sqrt(np.cos(np.radians(beta2BArr))) / (ZBarr ** 0.7))
sigmaMat = np.copy(etaMat)
b2Mat = np.copy(etaMat)
c2Mat = np.copy(etaMat)
c2mMat = np.copy(etaMat)
PrestMat = np.copy(etaMat)
VslipMat = np.copy(etaMat)
pressErrorMat = np.copy(etaMat)
MachExitMat = np.copy(etaMat)
WxMat = np.copy(etaMat)
Ctheta2Mat = np.copy(etaMat)
dh0SlipCorrMAt = np.copy(etaMat)
beta2flowMat= np.copy(etaMat)
# Declaring variables for off-design.py
D2 = 0
b2 = 0
U1t = 0
beta1 = 0






""" ------------------------------------ Beginning iteration ------------------------------------ """
print('\n')

dh0s = ( (settings.k * settings.R * settings.T00) / (settings.k - 1) ) * (settings.Pr **((settings.k - 1)/ settings.k) - 1)                     # Specific isentropic compression enthalpy, constant value [kJ/kg/K]
   
""" Calculating critical property of a rotating disk to use as constraint for RPM/rotational velocity/radiusettings. Disk will break at outermost point, therefore r2 and U2. """
omegaCritMultr1 = np.sqrt(3*settings.impellerTensileStrength/ settings.impellerDensity )/(1/settings.r1Divr2)
NcritMultr1 = omegaCritMultr1 *60 /(2*math.pi)
U2Crit = np.sqrt(2*settings.impellerTensileStrength/settings.impellerDensity)
U2 = U2Crit

""" Iterating to find a impeller tip velocity that satisfy material constraint.    """
while (U2) >= U2Crit and Ndes > 0:
    # U2omega
    # U1divr2 = 2*np.pi * settings.r1Divr2 * Ndes / 60     
    # secondOrderCoeff = [mu*(U2divr2**2), -np.mean(U1divr2)*Ctheta1, -Wx]
    # r2roots = np.roots(secondOrderCoeff)
    # r2 = max(r2roots)

    U1ti, W1ti, r1, Ctheta1, C1, T1, M1, P1, rho1, A1, U1t, Cm1, beta1, W1t, omega = inducerFlowWRTminimumRelativeVelocity(N=Ndes)
    U2 = (U1t/settings.r1Divr2)

    if (U2) > settings.bladeVelUpperLimit:
        break
    if (U2) >= U2Crit:
        Ndes -= 100
    r2 = U2/omega
    Ncrit = 60*U2Crit/(2*np.pi*r2)
    # U2 = U2divr2*r2


""" Iterating main functionality """
for iN in range(1): # To be iterated for rpm N

    for iz in range(len(ZBarr)):
        z = ZBarr[iz]
        debug=1         # breakpoint
        
        for ib in range(len(beta2BArr)):
            expCond = rhsExpLimitMat[ib, iz]
            sigma = sigmaWiesnerMat[ib, iz]           # slip factor
            etaStage = settings.etaStage0
            trueFalse2 = False
            np.rad2deg((beta2BArr[ib])[0])
            lookForBetterEta = 0
            

            D2 = 2*r2
            rh1 = settings.rhDivr1*r1
            NDivr1 = Ndes/r1


            trueFalse1 = False
            rhs = rhsExpLimitMat[ib, iz]
            
            if settings.r1Divr2 < rhs:
                trueFalse1 = True
            elif settings.r1Divr2 > rhs:
                trueFalse1 = True
                sigma = sigma* (1  -  (   (settings.r1Divr2-rhs)/(1-rhs)   )**3   )        # Correcting if ratio of radii is freater than epsilon limit
            else: 
                trueFalse1 = False
                continue

            trueFalseMat[ib, iz][0] = trueFalse1
            Vslip = (1 - sigma)*U2                              # slip velocity

            mu = sigma * settings.lambda2 / (settings.lambda2 - math.tan((beta2BArr[ib])[0]))                # Work input coefficient [-]    MSG: Can formula be adjusted? tan(-x=-tan(x)) 

            Ctheta2m = mu * U2                                # Absolute tangential exit velocity [m/s]         from work coefficient
            Cm2m = Ctheta2m / settings.lambda2                       # Absolute meridional exit velocity [m/s]         
            C2 = (Ctheta2m ** 2 + Cm2m ** 2) ** 0.5           # Absolute exit velocity [m/s]                    from pythagoras 

            # finding velocities related to slip
            # See slip velocity triangle in related document        
            Ctheta2ideal = U2 - Ctheta2m*np.tan((np.abs(( beta2BArr[ib] )[0])))   
            CTheta2Real = sigma*Ctheta2ideal         
            Vslip1 = Ctheta2ideal -CTheta2Real
            beta2flow = np.rad2deg( np.arctan( (Vslip + Cm2m*np.tan( np.abs(beta2BArr[ib])) )/Cm2m ) )

            


            dh0SlipCorrected =  U2*(CTheta2Real) - U1t*Ctheta1   # Alternative for finding fluid enthalpy change, for comparison

            while ( (trueFalseMat[ib, iz][1] == False) and (etaStage < settings.etaUpperLimit and etaStage > settings.etaLowerLimit) ):
                
                Wx = dh0s / etaStage                                     # Specific work [J/kg/K]
                workError = np.abs(Wx-dh0SlipCorrected)/Wx               # Comparing the two methods of finding enthalpy change
                wTest1 = Wx
                wTest2 = U2*Ctheta2m - U1t*Ctheta1

                Wx = dh0SlipCorrected
                

                # ------------------- Impeller outlet calculation. Stagnation properties denoted by zero. -------------------
                T02m = settings.T00 + Wx * (settings.k - 1) / (settings.k * settings.R)                                                     # Stagnation exit temperature [K]     , from dh=cp*Dt   (settings.k - 1) / (settings.k * settings.R) = 1/cp
                # P02m =( ( T02m/settings.T00 )**(settings.k/(settings.k-1)) )*settings.P00
                P02m = settings.P00 * ((etaStage * Wx * (settings.k - 1)  / (settings.k * settings.R * settings.T00)) + 1) ** (settings.k / (settings.k - 1))     # Exit Stagnation Pressure [Pa]       , from (... todo ...)


                T2m = T02m - (settings.k - 1) / (2 * settings.k * settings.R) * (C2 **2)                                             # Exit temperature [K]                            from stagnation temperature
                P2m = P02m / ((T02m / T2m) ** (settings.k / (settings.k - 1)))                                                # Exit pressure [Pa]                              from stagnation pressure, isentropic relation
                rho2m = P2m / (T2m * settings.R)                                                                       # Exit density [kg/m^3]                           from ideal gas law
                A2 = settings.mdot / (rho2m * Cm2m)                                                                    # Exit area [m^2]                                 from continuity
                b2 = A2 / (math.pi * D2)                                                                        # Depth of impeller exit [m]                      from geometry
                M2 = U2/np.sqrt(settings.k*settings.R*T02m)                                                                   # Impeller exit blade mach number 
                

                # ------------------- Diffuser Calculation -------------------
                P3 = P2m + settings.CpD * (P02m - P2m)               # Diffuser exit static pressure [Pa]
                C3 = C2 / settings.AR                                # Diffuser exit absolute velocity [m/s]
                P03 = P3 + 0.5 * rho2m * C3 ** 2              # Diffuser exit stagnation pressure [Pa]


                # ------------------- Overall performance -------------------
                etaiterate = ( (P03 / settings.P00) ** ((settings.k - 1) / settings.k) - 1 ) / ( (T02m / settings.T00) - 1 )        # Iterative stage efficiency [-], isothermal diffuser assumed?
                Prest = ((etaiterate * U2 ** 2 * mu) / (settings.Cp * T1) + 1) ** (settings.k / (settings.k - 1))          # Estimate of the pressure ratio, equation is validated
                PressureTestOuterLoop = (Prest-settings.Pr)/settings.Pr
                # print(etaiterate)

                """Checking if teration conditions are sustained"""
                if (etaiterate > settings.etaUpperLimit or etaiterate < settings.etaLowerLimit) or ( etaStage > settings.etaUpperLimit or etaStage < settings.etaLowerLimit ):
                    trueFalse2 = False
                elif abs(PressureTestOuterLoop) > settings.iterTol:
                    checkOuterLoop = False
                    trueFalse2 = False
                else:
                    trueFalse2 = True


                # trueFalse2 = True #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DISKUTER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                trueFalseMat[ib, iz][1] = trueFalse2

                # trueFalse2 = True #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DISKUTER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                if ib == 20 and iz ==20:
                    debug = 1

                if etaStage <= settings.etaLowerLimit or etaStage >= settings.etaUpperLimit:
                    break
                
                """ Updating values from nan's if iteration criterion is upheld"""
                if  ( trueFalse2 == True ) and U2 < settings.bladeVelUpperLimit:
                    # etaMat[ib, iz] = etaStage
                    etaMat[ib, iz] = etaiterate
                    pressErrorMat[ib, iz] = abs(PressureTestOuterLoop)
                    MachExitMat[ib, iz] = M2
                    b2Mat[ib, iz] = b2
                    c2Mat[ib, iz] = C2
                    c2mMat[ib, iz] = Cm2m
                    PrestMat[ib, iz] = Prest
                    VslipMat[ib, iz] = Vslip1
                    Ctheta2Mat[ib, iz] = Ctheta2m
                    WxMat[ib, iz] = Wx
                    dh0SlipCorrMAt[ib, iz] = dh0SlipCorrected
                    sigmaMat[ib, iz] = sigma
                    beta2flowMat[ib, iz] = beta2flow



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
betaMinErr = round(np.rad2deg(float(beta2BArr[IBmin][0])), 3)
zMinErr = round(float(ZBarr[IZmin]),3)
etaMinErr = round(float(etaMat[IBmin, IZmin]), 3)
effText = r"Minimum error using $Z_{B}$ = " + str(zMinErr) + r" and $\beta _{2b}$ = " + str(betaMinErr) + r" with efficiency $\eta$ = " +str(etaMinErr) 

# making nice text to go with plots
maxEta = np.nanmax(etaMat)
minEta = np.nanmin(etaMat)
countTrue = np.sum(np.all(trueFalseMat, axis=-1))
totalCases = len(ZBarr)*len(beta2BArr)


text1 = "\n" \
        r"Valid cases: " + str(countTrue) + '/'  +str(totalCases) + "\n" \
        r"Desired $PR^*$: " + str(settings.Pr) +  "\n" \
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
        r"$\frac{r_h}{r_1}{*}$ = "  +str(settings.rhDivr1) + r" $\frac{m}{m} $" +"    " \
        r"$\frac{r_1}{r_2}{*}$ = " + str(settings.r1Divr2) + r" $\frac{m}{m} $" +"\n" \
        r"$r_{t1}$ = " +str(round(r1, 3)) + r"$m$" +" \n" \
        r"$r_{h1}$ = " +str(round(rh1, 3)) + r"$m$" +" \n" \
        r"$r_{t2}$ = " +str(round(r2, 3)) + r"$m$" +" \n" \
        r"$N^{*}$ = " +str(round(Ndes, 1)) + r"$rpm$" +" \n" \
        r"$N_{crit}$ = " +str(round(Ncrit, 3)) + r"$rpm$" +" \n" \
        
text3 = "\n" \
        r"Proposed efficiency $\eta *=$  " + str(settings.etaStage0) + "\n" \
        r"Exit swirl number $\lambda_2 *$: " + str(settings.lambda20) +  "\n" \
        r"Tolerance*: " + str(round(settings.iterTol, 3)) + "   " \

""" Preparing input to plottting functions"""
systemVar = [settings.etaStage0,settings.lambda2, settings.iterTol, ZBarr, beta2BArr]
designParam = [r1, r2, rh1, Ncrit, Ndes]
flowVar = [settings.Pr, W1t, Cm1, U1t, U2, U2Crit]

Zflow = [Ctheta2Mat, pressErrorMat, etaMat, MachExitMat, PrestMat, trueFalseMat, WxMat, dh0SlipCorrMAt]
Zcompressor = [b2Mat, VslipMat, beta2flowMat]
Zsystem = [WxMat, dh0SlipCorrMAt,sigmaMat,etaMat, PrestMat]
Zvelocities = [c2Mat, Ctheta2Mat, c2mMat, MachExitMat, VslipMat]
text = [text1, text2, text3]

""" Plot of relative velocities """
fig, ax11 = plt.subplots()
plt.plot(settings.Cm1i, W1ti, '--g')
ax11.tick_params(axis='both', which='major', labelsize=10)
plt.xlabel(r"$C_{m1}  [m/s]$",fontsize=12)
plt.ylabel(r'$W_{1t} [m/s]$', fontsize=12)
# plt.title("Minimisation of W1t")
plt.grid()


plotCompressorParam(systemVar, Zcompressor, designParam, flowVar, text )
plotFlowConditions(systemVar, Zflow, designParam, flowVar, text)
plotVelocities(systemVar, Zvelocities, designParam, flowVar, text)
plotText(text)
plotSystemVariables(systemVar, Zsystem, designParam, flowVar, text)



# plt.show(block=False)
# show = logic.showGeometryPlots
plt.show()
