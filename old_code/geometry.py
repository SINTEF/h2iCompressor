"""
The following Python code takes a set of variables as inputs and calculates the geometry of an inducer, 
impeller and diffuser for a centrifugal compressor.

Reference: Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244
Author: Petter Resell (summer intern, 2024)

"""

### Import-------------------------------------------------------------------------------------------
import math
import numpy as np
import matplotlib.pyplot as plt
import settings


# Plotting files and functions:
from plot_scripts.plot_flow import plotFlowConditions
from plot_scripts.plot_system import plotSystemVariables
from plot_scripts.plot_compressor import plotCompressorParam
from plot_scripts.plot_velocities import plotVelocities
from plot_scripts.plot_text import plotText
# from find_flow import inducerFlowWRTminimumRelativeVelocity
from pressure_test import pressureOverUnderEstimate
import geometry_system_functions

plt.rcParams.update(plt.rcParamsDefault)            # For plotting
plt.rcParams.update({'font.size': 15})

print('\n' + 'Running Geometry.py')

Ndes = settings.N0                                   # Rotational speed [rpm]     
etaStage = settings.etaStage0                        # Isentropic Stage efficiency [-]           


### ------------------- Initializing arrays and matrices for iteration and later avaluation -------------------

""" ZB is row vector, varies along row, constant at column. beta2B is column vector, varies along column height, constant along row distance """
ZBarr = np.arange(settings.bladeMin, settings.bladeMax + 1, 1)                     # Blade number          
beta2BArr = np.radians(np.arange( settings.beta2Bmax, settings.beta2Bmin, 1))      # discharge angle
beta2BArr = beta2BArr[:, np.newaxis]                                                                 # Flipping from row to column vector to make matrix on next lines

""" Making matrices for iteration later. Filling with nans that are only replaced if all conditions are met. """
rhsExpLimitMat = np.exp(-8.16* np.cos( beta2BArr)/ ZBarr)                                                                   # Matrix for epsilon_limit from wiesner condition           
etaMat = np.array([[ np.nan for _ in range(np.shape(rhsExpLimitMat)[1])] for _ in range(np.shape(rhsExpLimitMat)[0])] )     # Matrix for efficiency. Replace by fill matrix with same shape as rhsExpLimitMAt with nans 
sigmaWiesnerMat = 1 - (np.sqrt(np.cos(np.radians(beta2BArr))) / (ZBarr ** 0.7))                                             # Matrix for wiesner slip factor 
sigmaMat = np.copy(etaMat)                                                                                                  # Matrix for slip factor found to be valid
b2Mat = np.copy(etaMat)                                                                                                     # Matrix for impeller exit cylinder height
c2Mat = np.copy(etaMat)                                                                                                     # Matrix for impeller absolute discharge velocity
c2mMat = np.copy(etaMat)                                                                                                    # Matrix for meridonal component of impeller discharge velocity
PrestMat = np.copy(etaMat)                                                                                                  # Matrix for pressure estimate
VslipMat = np.copy(etaMat)                                                                                                  # Matrix for slip velocity
pressErrorMat = np.copy(etaMat)                                                                                             # Matrix for pressure error
MachExitMat = np.copy(etaMat)                                                                                               # Matrix for impeller mach number
WxMat = np.copy(etaMat)                                                                                                     # Matrix for compression work
Ctheta2Mat = np.copy(etaMat)                                                                                                # Matrix for angular component of discharge velocity
dh0SlipCorrMAt = np.copy(etaMat)                                                                                            # Matrix for work found by slip corrected euler equation
beta2flowMat= np.copy(etaMat)                                                                                               # Matrix for discharge flow angle found by slip relations

# Declaring variables for off-design.py
U1t = 0
beta1 = 0


### Inducer calculation----------------------------------------------------------------------------------
"""
For each value of Cm1 in settings.Cm1i, we calculate the corresponding values of C1, T1, M1, P1, rho1, A1, r1, U1t, W1t, Cm1 and beta1.
    The value of Cm1 will affect the values of U1t and Ctheta1 (through ex. r1), which in turn affects the value of W1t. Isentropic relations 
        linking the static to the total conditions. Velocity triangles are also applied. 
            By plotting W1t against Cm1, we then find the value of Cm1 that minimises W1t.

"""
Ctheta1i = settings.Cm1i * math.tan(math.radians(settings.alpha1))                                    # Inducer absolute tangential velocity [degrees]          Trigonometry
C1i = (Ctheta1i ** 2 + settings.Cm1i ** 2) ** 0.5                                                              # Inducer Absolute velocity C1 [m/s]                      Pythagoras theorem from velocity triangle 
T1i = settings.T00i - (C1i ** 2) / (2 * settings.Cp)                                                  # Inducer temperature [K]                                 Stagnation temperature relation             
M1i = C1i / ((settings.k * settings.R * T1i) ** 0.5)                                                  # Inducer Mach number [-]                                 by definition
P1i = settings.P00 * (T1i / settings.T00) ** (settings.k / (settings.k - 1))        # Inducer pressure [Pa]                                   Isentropic relation                        
rho1i = P1i / (settings.R * T1i)                                                                               # Inducer density  [kg/m^3]                               Assuming ideal gas                          
A1i = settings.mdot / (rho1i * settings.Cm1i * (1 - settings.B1))                            # Inducer flow area [m^2]                                 Continuity                                                        
rt1i = (A1i / (math.pi* (1 - settings.rhDivr1**2) )) ** 0.5                                                    # Inducer tip radius [m]                                  Geometry


""" This function takes a rotational speed N and findes the relative velocity. This is then minimized to increase efficiency. 
        Inducer properties and velocities are then found at the point where the relative velocity is the smallest. 
            This is at w1ti_index_min """
def inducerFlowWRTminimumRelativeVelocity(N):
    U1ti = 2 * math.pi * rt1i * N / 60                                           # Inducer blade tip speed [m/s]
    W1ti = (settings.Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5         # Inducer relative velocity [m/s]
    w1ti_index_min = np.argmin(W1ti)                                                 # Index of smallest relative velocity [m/s]
    r1 = rt1i[w1ti_index_min]                                                        # Inducer tip radius [m]
    Ctheta1 = Ctheta1i[w1ti_index_min]                                               # Inducer angular velocity [m/s]
    C1 = C1i[w1ti_index_min]                                                         # Inducer velocity [m/s]
    T1 = T1i[w1ti_index_min]                                                         # Inducer temperature [K]
    M1 = M1i[w1ti_index_min]                                                         # Inducer Mach number [-]
    P1 = P1i[w1ti_index_min]                                                         # Inducer pressure [Pa]
    rho1 = rho1i[w1ti_index_min]                                                     # Inducer density [kg/m3]
    A1 = A1i[w1ti_index_min]                                                         # Inducer flow area [m3]
    U1t = U1ti[w1ti_index_min]                                                       # Inducer tip speed [m/s]
    Cm1 = settings.Cm1i[w1ti_index_min]                                     # Inducer meridonal velocity [m/s]
    beta1 = math.degrees(math.atan((U1t - Ctheta1) / Cm1))                       # Inducer relative velocity angle [deg]
    W1t = W1ti[w1ti_index_min]                                                       # Inducer relative velocity [m/s]
    omega = U1t/r1                                                               # Angular velocity [rad/s]          =(2*np.pi*N/60)
    return U1ti,W1ti,r1,Ctheta1,C1,T1,M1,P1,rho1,A1,U1t,Cm1,beta1,W1t,omega


### ------------------------------------  Impeller calculation  --------------------------------------

""" The isentropic enthalpy demand is used to find the approximate required work. This is done through iterating the isentropic efficiency.
    Furthermore the  Work is used to find the outlet stagnation temperature which in turn is used to find the outlet stagnation pressure 
    from isentropic relations. Furthermore, geometrical relations are applied. """

""" If the pressure estimate Prest initially over-estimates the pressure ratio by more than 5%, then the the error will only grow for increasing efficiency
        since the pressure ratio increase for increasing efficiency. Hence its necessary to check if the pressure estimate is an 
            under- or over-estimate before adjusting the efficiency eta. From then the efficiency can be either increased or decreased. """


""" Specific isentropic compression enthalpy, constant value [kJ/kg/K]. dh0s is constant throughout the full iteration scheme. """
dh0s = ( (settings.k * settings.R * settings.T00) / (settings.k - 1) ) * (settings.Pr **((settings.k - 1)/ settings.k) - 1)                     



""" Calculating critical property of a rotating disk to use as constraint for RPM/rotational velocity/radiusettings. Disk will break at outermost point, therefore r2 and U2. """                                                                                                 
U2Crit = np.sqrt(2*settings.impellerTensileStrength/settings.impellerDensity)      # Applying tensile strength of disk. Titanium used.        
U2 = U2Crit        
                                                                                  # initializing for loop 

""" Iterating to find a impeller tip velocity that satisfy material constraint. Taking the rotational speed set in settings.py and
         checking if the corresponding rotational velocity U2 is bigger than the critical rotational velocity. If it is then 
            inducerFlowWRTminimumRelativeVelocity(N=Ndes) finds all inducer flow properties for the new N and U2. For the first iteration U2=U_crit 
              so it runs at least one iteration anyways. If the new U2 is found to be to large then it is lowered by increments of 1000.
                This could be replaced by seting the rotational velocity to for example 70% or 80% of the critical one.  """
while (U2) >= U2Crit and Ndes > 0:

    U1ti, W1ti, r1, Ctheta1, C1, T1, M1, P1, rho1, A1, U1t, Cm1, beta1, W1t, omega = inducerFlowWRTminimumRelativeVelocity(N=Ndes)
    U2 = (U1t/settings.r1Divr2)

    if (U2Crit) > settings.bladeVelUpperLimit:      # Can skip this break to avoid crashing code
        break
    if (U2) >= U2Crit:
        Ndes -= 1000


""" Can find design parameters after finding U2"""
r2 = U2/omega                                           # Impeller tip radius
D2 = 2*r2                                               # Impeller tip diameter
rh1 = settings.rhDivr1*r1                      # Hub radius
Ncrit = 60*U2Crit/(2*np.pi*r2)                          # Critical rpm for the newly found impeller tip radius. 


""" Function to check if condition for the use of wiesner slip factor is upheld. If its not upheld then slip factor 
        is updated by the correction equation. If the condition is upheld it just returns the same slip factor. 
         The Function is implemented when iterating through all blade numbers and blade angles below.  """
def checkUpdateSlipFactorSigma(epsilonLimit, slipFactor):
    if settings.r1Divr2 >= epsilonLimit:
        return slipFactor* (1  -  (   (settings.r1Divr2-epsilonLimit)/(1-epsilonLimit)   )**3   )
    else:
        return slipFactor


""" Iterating main functionality of design procedure. Looping through each blade number and blade angle to 
        capture all combinations.  """    

for iz in range(len(ZBarr)):
    # debug=1         # breakpoint for debugging
    
    for ib in range(len(beta2BArr)):

        sigma = sigmaWiesnerMat[ib, iz]                 # slip factor for given blade number and blade angle
        etaStage = settings.etaStage0          # resetting the guessed efficiency for each new combination of blade number and blade angle
        trueFalseCheck = False                          # iteration control for while loop
        

        """ Handling slip factor. Updating if the conditions for wiesner relation is not upheld. """
        epsilonLimit = rhsExpLimitMat[ib, iz]                                                          # right hand side of wiesner slip factor condition
        sigma = checkUpdateSlipFactorSigma(epsilonLimit = epsilonLimit, slipFactor = sigma)            # updating slip factor wrt wiesner condition
        
        """ impellerOutletVelocities() finds the impeller outlet velocities and work output coeff. """
        workInputCoeff, Ctheta2m, Cm2m, C2 = geometry_system_functions.impellerOutletVelocities(slipFactor = sigma, beta2B = (beta2BArr[ib])[0], U2 = U2)          # Finding impeller outlet velocities


        """ The following five lines are made to find porperties of the slip for a given blade number and blade angle.
                None of the resulting variables are applied in further calculations, but could easily be. They are included
                 for further development and demonstration of slip.  """
        Ctheta2ideal = U2 - Ctheta2m*np.tan((np.abs(( beta2BArr[ib] )[0])))   
        CTheta2Real = sigma*Ctheta2ideal         
        Cslip1 = Ctheta2ideal -CTheta2Real
        beta2flow = np.rad2deg( np.arctan( (Cslip1 + Cm2m*np.tan( np.abs(beta2BArr[ib])) )/Cm2m ) )
        dh0SlipCorrected =  U2*(CTheta2Real) - U1t*Ctheta1   # Alternative for finding fluid enthalpy change, for comparison




        """ The while loop iterates for a given combination of blade number and blade angle. It iterates the efficiency demand
                while the estimated pressure ratio is not within the given tolerance of the desired pressure ratio. It is updated in 
                 the function called pressureOverUnderEstimate. The guess efficiency is increased or decreased varying if the 
                  oressure ratio is under- or overestimated. """
        while ( (trueFalseCheck == False) and (etaStage < settings.etaUpperLimit and etaStage > settings.etaLowerLimit) ):
            
            # ------------------- Updating work with the new isentropic efficiency -------------------
            Wx = dh0s / etaStage                                     # Specific work [J/kg/K]
            workError = np.abs(Wx-dh0SlipCorrected)/Wx               # Comparing the two methods of finding enthalpy change
            # Wx = dh0SlipCorrected                                   # Work set to be given by slip corrected euler equation 
            

            """ Impeller outlet calculation. Stagnation properties denoted by zero. Isentropic relations applied on next four rows. """
            T02m = settings.T00 + Wx * (settings.k - 1) / (settings.k * settings.R)                                                     # Stagnation exit temperature [K]     , from dh=cp*Dt   
            P02m = settings.P00 * ((dh0s * (settings.k - 1)  / (settings.k * settings.R * settings.T00)) + 1) ** (settings.k / (settings.k - 1))     # Exit Stagnation Pressure, isentropic [Pa]      
            T2m = T02m - (settings.k - 1) / (2 * settings.k * settings.R) * (C2 **2)             # Exit temperature [K]                            from stagnation temperature
            P2m = P02m / ((T02m / T2m) ** (settings.k / (settings.k - 1)))                                # Exit pressure [Pa]                              from stagnation pressure, isentropic relation
            rho2m = P2m / (T2m * settings.R)                                                                       # Exit density [kg/m^3]                           from ideal gas law
            A2 = settings.mdot / (rho2m * Cm2m)                                                                    # Exit area [m^2]                                 from continuity
            b2 = A2 / (math.pi * D2)                                                                                        # Impeller exit cylinder height [m]               from geometry
            M2 = U2/np.sqrt(settings.k*settings.R*T02m)                                                   # Impeller exit blade mach number                 from definition
            

            # ------------------- Finding diffuser properties -------------------
            P3, P03, C3 = geometry_system_functions.diffuserFlow(P2=P2m, P02=P02m, rho2=rho2m, C2=C2)


            # ------------------- Overall performance -------------------
            etaiterate, Prest, PressureTestOuterLoop = geometry_system_functions.systemTotalPerformance(P03, T02m, U2, T1, workInputCoeff)
            
            """ Checking if iteration conditions are sustained. the iterated efficiency is only accepted if
                    the pressure estimate error is less than the iteration tolerance and if the iterated
                     efficiency is within the reasonable limits etaUpperLimit and etaLowerLimit """
            if (etaiterate > settings.etaUpperLimit or etaiterate < settings.etaLowerLimit) or ( etaStage > settings.etaUpperLimit or etaStage < settings.etaLowerLimit ):
                trueFalseCheck = False
            elif abs(PressureTestOuterLoop) > settings.iterTol:
                checkOuterLoop = False
                trueFalseCheck = False
            else:
                trueFalseCheck = True

            # trueFalseCheck = True


            if etaStage <= settings.etaLowerLimit or etaStage >= settings.etaUpperLimit:
                break
            
            """ Updating values from nan's if iteration criterion is upheld"""
            if  ( trueFalseCheck == True ) and U2 < settings.bladeVelUpperLimit:
                etaMat[ib, iz] = etaiterate
                pressErrorMat[ib, iz] = abs(PressureTestOuterLoop)
                MachExitMat[ib, iz] = M2
                b2Mat[ib, iz] = b2
                PrestMat[ib, iz] = Prest
                WxMat[ib, iz] = Wx
                dh0SlipCorrMAt[ib, iz] = dh0SlipCorrected
                c2Mat[ib, iz] = C2
                c2mMat[ib, iz] = Cm2m
                Ctheta2Mat[ib, iz] = Ctheta2m
                VslipMat[ib, iz] = Cslip1
                sigmaMat[ib, iz] = sigma
                beta2flowMat[ib, iz] = beta2flow



            """ Updating efficiency for next iteration"""
            etaStage = pressureOverUnderEstimate(PressureTestOuterLoop, etaStage)
    
    
        #     debug = 1       # Breakpoint for debugging
        # debug =1            # Breakpoint for debugging
            



# ------------- Preparing data for nice plotting ------------ 
""" Plotting functions plotCompressor.py, plotFlow.py, plotSystem.py, plotVelocities.py and plotText.py are 
        called from logic.py but prepared here since data is produced here."""


""" making nice text to go with plots """
maxEta = np.nanmax(etaMat)
minEta = np.nanmin(etaMat)
countTrue = np.count_nonzero(~np.isnan(etaMat))
totalCases = len(ZBarr)*len(beta2BArr)

text1 = "\n" \
        r"Valid cases: " + str(countTrue) + '/'  +str(totalCases) + "\n" \
        r"$\eta _{max}$ = " + str(round(maxEta, 3)) + "   " \
        r"$\eta _{min}$ = " + str(round(minEta, 3)) + "\n" \
        r"$W_{1t}= $" + str(round(W1t, 3)) + "m/s   " \
        r"$C_{m1}= $" + str(round(Cm1, 3)) + "m/s"+"\n" \
        r"$U_{1t}= $" + str(round(U1t, 3)) + "m/s"+"   " \
        r"$M_{1f}= $" + str(round(M1, 3)) + "\n" \
        r"$U_{2t}$ = "  +str(round(U2, 3)) + "m/s"+"   " \
        r"$U_{2t,crit}= $" + str(round(U2Crit, 3)) +"m/s"+ "\n" \

text2 = "\n" \
        r"$r_{t1}$ = " +str(round(r1, 3)) + r"$m$" +" \n" \
        r"$r_{h1}$ = " +str(round(rh1, 3)) + r"$m$" +" \n" \
        r"$r_{t2}$ = " +str(round(r2, 3)) + r"$m$" +" \n" \
        r"$N_{crit}$ = " +str(round(Ncrit, 3)) + r"$rpm$" +" \n" \
        
text3 = "\n" \
        r"Desired $PR^*$: " + str(settings.Pr) +  "\n" \
        r"$\frac{r_h}{r_1}{*}$ = "  +str(settings.rhDivr1) + r" $\frac{m}{m} $" +"    " \
        r"$\frac{r_1}{r_2}{*}$ = " + str(settings.r1Divr2) + r" $\frac{m}{m} $" +"\n" \
        r"$N^{*}$ = " +str(round(Ndes, 1)) + r"$rpm$" +" \n" \
        r"Proposed efficiency $\eta *=$  " + str(settings.etaStage0) + "\n" \
        r"Exit swirl number $\lambda_2 *$: " + str(settings.lambda20) +  "\n" \
        r"Tolerance*: " + str(round(settings.iterTol, 3)) + "   " \

""" Preparing input to plottting functions"""
""" Common input for all plotting functions"""
systemVar = [settings.etaStage0,settings.lambda2, settings.iterTol, ZBarr, beta2BArr]
designParam = [r1, r2, rh1, Ncrit, Ndes]
flowVar = [settings.Pr, W1t, Cm1, U1t, U2, U2Crit]

""" One Z(...)-array goes to one plot function. For example: Zflow goes to plotFlow.py, Zsystem goes to plotSystem.py etc."""
Zflow = [Ctheta2Mat, pressErrorMat, etaMat, MachExitMat, PrestMat, WxMat, dh0SlipCorrMAt]     # Goes to plotFlow.py
Zcompressor = [b2Mat, VslipMat, beta2flowMat]                                                               # Goes to plotCompressor.py
Zsystem = [WxMat, dh0SlipCorrMAt,sigmaMat,etaMat, PrestMat]                                                 # Goes to plotSystem.py
Zvelocities = [c2Mat, Ctheta2Mat, c2mMat, MachExitMat, VslipMat]                                            # Goes to plotVelocities.py
text = [text1, text2, text3]                                                                                # Goes to plotText.py

""" Plot of relative velocities """
# fig, ax11 = plt.subplots()
# plt.plot(settings.Cm1i, W1ti, '--g')
# ax11.tick_params(axis='both', which='major', labelsize=10)
# plt.xlabel(r"$C_{m1}  [m/s]$",fontsize=12)
# plt.ylabel(r'$W_{1t} [m/s]$', fontsize=12)
# # plt.title("Minimisation of W1t")
# plt.grid()

print('Critical RPM: ' + str(round(Ncrit, 2)) )
print('Applied RPM: ' + str(round(Ndes, 2)) )
print('Rotational velocity: ' + str(round(U2, 2)))
print('rh: ' + str(round(rh1, 3)) + 'm' )
print('r1: ' + str(round(r1, 3)) + 'm')
print('r2: ' + str(round(r2, 3)) + 'm' )

if countTrue ==0:
    raise Exception(f"combination of rh/r1, r1/r2 and N gave zero valid cases. RPM to low to achieve desired pressure ratio!")

print('Geometry.py successfully run. \n')