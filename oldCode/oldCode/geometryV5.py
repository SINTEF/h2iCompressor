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
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams.update({'font.size': 15})

### Variables-------------------------------------------------------------------------------------------
# Inlet flow parameters
mdot = 30                         # Mass flow rate [kg/s]   , mdot = ( (300* 10**3 )/( 60**2 ) )/6  

Ndes = 120000                           # Rotational speed [rpm]
Ndes = 60000                          # Rotational speed [rpm]
Ndes = 10000                          # Rotational speed [rpm]


Cm1i = np.arange(10, 600.5, 0.5)       # Inlet Absolute meridional velocity [m/s]

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

# Steel properties
steelTensileStrength = 500* 10**6       # Approximate ultimate tensile strength UTS [Pa]
steelDensity = 8000                     # [kg/m3]

# Impeller exit parameters
lambda2 = 2                 # Exit swirl parameter                                          , check out
lambda20 = lambda2          # Exit swirl parameter, used for iteration
etaStage = 0.6             # Isentropic Stage efficiency [-]           
etaStage0 = etaStage                               
etaiterate = etaStage       # Isentropic stage efficiency used for iteration [-]

# Generall compressor parameters
# inlet conditions
P00 = 30 * (10**2) * (10**3)                    # Inlet stagnation pressure [Pa]
T00 = 293                                       # Inlet stagnation temperature [K]
Pr = 1.24                   # Pressure ratio [-]

alpha1 = 0                                      # Absolute inlet velocity angle [degrees]
rh1 = 0.2                                      # Hub radius [m]
B1 = 0.04                                       # Boundary layer blockage [-],        substance viscocity dependant
B2 = 0.04                                       # Boundary layer blockage [-],        substance viscocity dependant NEED GREATER BLOCKAGE FOR OUTLET
AR = 2.5                                        # inlet/outlet area ratio of the diffuser [-]         
T00i = np.full(len(Cm1i), 293)                  # Inlet stagnation temperature [K]


### Inducer calculation----------------------------------------------------------------------------------
"""
For each value of Cm1 in Cm1i, we calculate the corresponding values of C1, T1, M1, P1, rho1, A1, r1, U1t, W1t, Cm1(?) and beta1.
The value of Cm1 will affect the values of U1t and Ctheta1 (through ex. r1), which in turn affects the value of W1t.
By plotting W1t against Cm1, we then find the value of Cm1 that minimises W1t.

r1 inlet tip radius ???

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
rhDivr1=0.5                                                    # Ratio commonly given
rt1i = (A1i / (math.pi* (1 - rhDivr1**2) )) ** 0.5 



# def findRelativeInletVelocity(Ndes):
U1ti = 2 * math.pi * rt1i * Ndes / 60                       # Inlet blade tip speed [m/s]                    eqn (60) rapport
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

def pressureOverUnderEstimate( outerTest, etaStage):
    if outerTest == np.nan:
        debug = 0                   # breakpoint for debugging
    if  outerTest > 0:
        etaST = etaStage - 0.002                  # todo: need better way of iterating etaStage
    elif  outerTest < 0:
        etaST = etaStage + 0.002

    return etaST

### -----------------------  Iteration control initialization -----------------------------------
etaLowerLimit = 0.25                        # Lowest efficiency allowed
etaUpperLimit = 0.95
#  etaUpperLimit = 1                          # Highest efficiency allowed
bladeVelUpperLimit = 1200                   # Highest blade velocity allowed
bladeVelLowerLimit = 0                      # lowest blade velocity allowed
betamax = -45                               # "Maximum" beta iterated over
betamin = 0                                 # "Minimum" beta iterated over  
bladeMin = 1
bladeMax = 30
iterTol = 0.015                             # Outer loop condition tolerance

lambda2Arr = np.arange(0.2)
ZBarr = np.arange(1, 31, 1)                                                             
beta2bArr = np.arange(np.radians(betamax), np.radians(betamin), np.radians(1))
beta2bArr = beta2bArr[:, np.newaxis]
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


pressErrorMat = np.copy(etaMat)
Narr = np.arange(5000, 100000, 1000)
for iz in range(len(ZBarr)):
    for ib in range(len(beta2bArr)):
        beta2ZBmat[ib, iz] = (ZBarr[iz], beta2bArr[ib][0])

NcritArr = np.copy(etaMat)
NArr = np.copy(etaMat)
MachExitArr = np.copy(etaMat)

breakCheck = False


r1Dr2 = np.arange(0.8, 0.3, -0.1)       # r1 divided by r2 ,making r2 increase through array progression
U1divr2 = math.pi * r1Dr2 * Ndes / 30         

omegaCritMultr1 = np.sqrt(3*steelTensileStrength/ steelDensity )/(1/r1Dr2)
NcritMultr1 = omegaCritMultr1 *60 /(2*math.pi)
U2Crit = np.sqrt(2*steelTensileStrength/steelDensity)
U2 = U2Crit

# Narr = np.zeros( (len(Ncrit), int(stepsN) ) )

# for iN in range(len(Ncrit)):
#     matrixRow = np.linspace(1000, Ncrit[iN], stepsN)
#     Narr[iN, :len(matrixRow)] = matrixRow

# NarrIterate = np.arange(1000, Ncrit)

""" ------------------------------------ Beginning iteration ------------------------------------ """
print('\n')

countcrit=0

dh0s = ( (k * R * T00) / (k - 1) ) * (Pr **((k - 1)/ k) - 1)                     # Specific isentropic compression enthalpy, constant value [kJ/kg/K]
   

for iN in range(1): # To be iterated for rpm N

    for iz in range(len(ZBarr)):
        z = ZBarr[iz]
        debug=1         # breakpoint
        
        for ib in range(len(beta2bArr)):
            expCond = rhsExpLimitMat[ib, iz]
            sigma = sigmaWmat[ib, iz]           # slip factor
            etaStage = etaStage0
            etaiterate = etaStage0
            trueFalse1 = False
            trueFalse2 = False
            b = np.rad2deg((beta2bArr[ib])[0])
            lookForBetterEta = 0
            # while ( (trueFalseMat[ib, iz][1] == False) and (etaiterate < etaUpperLimit and etaiterate > etaLowerLimit) ):
            while ( (trueFalseMat[ib, iz][1] == False) and (etaStage < etaUpperLimit and etaStage > etaLowerLimit) ):
                mu = sigma * lambda2 / (lambda2 - math.tan((beta2bArr[ib])[0]))                # Work input coefficient [-]    MSG: Can formula be adjusted? tan(-x=-tan(x)) 
                # Todo: Mener etaStage i uttrykket over skal bort burde være:
                # P02m = P00*(T02m/T00)**(k/(k-1))

                Wx = dh0s / etaStage                                                             # Specific work [J/kg/K]
                T02m = T00 + Wx * (k - 1) / (k * R)                                              # Stagnation exit temperature [K]     , from dh=cp*Dt   (k - 1) / (k * R) = 1/cp
                P02m = P00 * ((etaStage * Wx * (k - 1) / (k * R * T00)) + 1) ** (k / (k - 1))    # Exit Stagnation Pressure [Pa]       , from (... todo ...)

                # U2 = ((U1t * Ctheta1 + Wx) / mu) ** 0.5    # Exit Blade Speed [m/s]
                # D2 = 60 * U2 / (math.pi * Ndes)              # Exit Diameter [m]
                # r2 = 0.5*60 * U2 / (math.pi * Ndes)              # Exit Diameter [m]
                
                U2divr2 = np.pi *Ndes /30
                U1divr2 = math.pi * r1Dr2 * Ndes / 30     
                secondOrderCoeff = [mu*(U2divr2**2), -np.mean(U1divr2)*Ctheta1, -Wx]
                r2roots = np.roots(secondOrderCoeff)
                r2 = max(r2roots)
                U2 = U2divr2*r2
                U1 = np.mean(U1divr2)*r2

                r1 = 30*U1/(np.pi*Ndes)         # r2* r1Dr2
                D2 = 2*r2
                rh1 = np.sqrt(r1**2 - (A1/np.pi))
                NDivr1 = Ndes/r1

  

                # correcting U2
                omegaCrit = np.sqrt(3*steelTensileStrength/(steelDensity * (r2)**2) )
                omegaCrit = (np.sqrt(2*steelTensileStrength/( steelDensity * r2**2 )) )
                Ncrit = omegaCrit*60/(2*math.pi)


                Vslip = (1 - sigma)/U2
                
                r1tDivr2t = r1/(r2)
                r1r2Mat[ib, iz] = r1tDivr2t

                trueFalse1 = False
                rhs = rhsExpLimitMat[ib, iz]
                if (etaiterate > 1 and etaiterate < 0):
                    trueFalse1 = False
                    break
                elif r1tDivr2t < rhs:
                    trueFalse1 = True
                else: 
                    trueFalse1 = False
                    break

                trueFalseMat[ib, iz][0] = trueFalse1

                #iterer på eta og finn max eta som oppfyller begge

                # U2 = D2 * math.pi * Ndes / 60
                    # !Dont touch!
                Ctheta2m = mu * U2                                # Absolute tangential exit velocity [m/s]         from work coefficient
                Cm2m = Ctheta2m / lambda2                         # Absolute meridional exit velocity [m/s]         

                muTest = Cm2m/U2
                phi2 = Cm2m/U2                                    # Flow coeff. [-]
                # Cm2m = (U2 - Ctheta2m)/math.tan(beta2b)         # Absolute meridional exit velocity [m/s] alternative 
                # sigma = 1 - (math.sqrt(math.cos(math.radians(beta2b))) / ( (ZB ** 0.7)*(1 - phi2 * math.tan(beta2b)) ))  # Alternative for slip factor iterere til kopnvergens????

                C2 = (Ctheta2m ** 2 + Cm2m ** 2) ** 0.5           # Absolute exit velocity [m/s]                    from pythagoras 
                T2m = T02m - (k - 1) / (2 * k * R) * (C2 **2)     # Exit temperature [K]                            from stagnation temperature
                P2m = P02m / ((T02m / T2m) ** (k / (k - 1)))      # Exit pressure [Pa]                              from stagnation pressure
                rho2m = P2m / (T2m * R)                           # Exit density [kg/m^3]                           from ideal gas law
                A2 = mdot / (rho2m * Cm2m)                        # Exit area [m^2]                                 from continuity
                b2 = A2 / (math.pi * D2)                          # Depth of impeller exit [m]                      from geometry
                M2 = U2/np.sqrt(k*R*T02m)
                

                # Diffuser Calculation----------------------------------------------------------------------------------
                etad = 0.85                                 # Estimated diffuser efficiency
                CpDi = 1 - 1 / (AR ** 2)                    # Ideal pressure recovery coefficient
                CpD = etad * CpDi                           # Pressure recovery coefficient
                P3 = P2m + CpD * (P02m - P2m)               # Diffuser exit static pressure [Pa]
                C3 = C2 / AR                                # Diffuser exit absolute velocity [m/s]
                P03 = P3 + 0.5 * rho2m * C3 ** 2            # Diffuser exit stagnation pressure [Pa]


                ### Overall performance----------------------------------------------------------------------------------
                """etaiterate becomes very small since T02m becomes so large. """
                etaiterate = ( (P03 / P00) ** ((k - 1) / k) - 1 ) / ( (T02m / T00) - 1 )        # Iterative stage efficiency [-], isothermal diffuser assumed?
                Prest = ((etaiterate * U2 ** 2 * mu) / (Cp * T1) + 1) ** (k / (k - 1))      # Estimate of the pressure ratio, equation is validated


                
                PressureTestOuterLoop = (Prest-Pr)/Pr
                if PressureTestOuterLoop > 0:
                    debug = 0
                    
                if (etaiterate > etaUpperLimit or etaiterate < etaLowerLimit) or ( etaStage > etaUpperLimit or etaStage < etaLowerLimit ):
                    trueFalse2 = False
                elif abs(PressureTestOuterLoop) > iterTol:
                    checkOuterLoop = False
                    trueFalse2 = False
                else:
                    trueFalse2 = True


                check = trueFalseMat[ib, iz][1]
                trueFalseMat[ib, iz][1] = trueFalse2

                if etaStage <= etaLowerLimit or etaStage >= etaUpperLimit:
                    break
                
                if ( trueFalse1 == True ) and ( trueFalse2 == True ) and U2 < bladeVelUpperLimit:
                    etaMat[ib, iz] = etaStage
                    U2Mat[ib, iz] = U2
                    r2Mat[ib, iz] = r2
                    pressErrorMat[ib, iz] = abs(PressureTestOuterLoop)
                    NcritArr[ib, iz] = Ncrit
                    MachExitArr[ib, iz] = M2
                    rh1Mat[ib, iz] = rh1
                    r1Mat[ib, iz] = r1
                    h2Mat[ib, iz] = b2
                    r2Mat[ib, iz] = r2
                    c2Mat[ib, iz] = C2
                    lookForBetterEta +=1
                        
                    if NDivr1 < np.mean(NcritMultr1):
                        countcrit +=1

                etaStage = pressureOverUnderEstimate( PressureTestOuterLoop, etaStage)
        
        
                debug = 1       #breakpoint
            debug =1
            


            


"""If the pressure estimate Prest initially over-estimates the pressure ratio by more than 5%, then the the error will only grow for increasing efficiency
    since the pressure ratio increase for increasing efficiency. 6
    Hence its neeeded to check if the pressure estimate is an under- or over-estimate before adjusting the efficiency eta. """







#-------------------- Code Outputs ------------------------------------------------------------
digits = 5
"""
print('--------------------------- Code Outputs ---------------------------')
print("Pressure ratio =", round(P03 / P00, digits) )
print("Pressure ratio error =", round(Pr - P03 / P00, digits))
print("Efficiency =", round(etaiterate, digits))
print("Efficiency error =", round( (abs(etaiterate - etaStage)), 5) )
print("Work =", round(Wx * mdot, digits-3), "J/kg")
print("Mass flow rate =", round(mdot, digits), "kg/s")
print("Exducer diameter =", round(D2 * 1000, digits), "mm")
print("Inducer Diameter =", round(r1 * 2000, digits), "mm")
print("Exit blade speed =", round(U2, digits))
print("Number of iterations of outer loop: ", int(countOuterLoop-1))


# if r1 / (0.5 * D2) < math.exp(- 8.16 * math.cos(beta2b) / ZB):     # Check for sigma validity
#     print("sigma valid")
# else:
#     print("sigma invalid!")

data1 = [round(P03 / P00, digits), round(Pr - P03 / P00, digits), round(etaiterate, digits), round( (abs(etaiterate - etaStage)), digits), 
          round(Wx * mdot, digits-3), round(mdot, digits), round(D2 * 1000, digits),round(r1 * 2000, digits),int(round(U2, digits)), P00* 10**-5, P03* 10**-5 ]
col1 = ["Pressure ratio", "Pressure ratio error", "Efficiency", "Efficiency error", "Work", "Mass flow rate","Exducer diameter", 
        "Inducer Diameter", "Exit blade speed", "Inlet Pressure", "Outlet Pressure"]
units = ['-', '-', '-', '-', 'Watt', 'kg/s', 'mm', 'mm', '?', '[bar]', '[bar]']
data = {"*Variable*    ":col1, "      *Value*":data1, "      *Unit*":units}
df = pd.DataFrame(data)
table = df.to_string(index=False)
print(table)
print('\n')
"""


maxEta = np.nanmax(etaMat)
minEta = np.nanmin(etaMat)
countTrue = np.sum(np.all(trueFalseMat, axis=-1))
totalCases = len(ZBarr)*len(beta2bArr)
print('Maximum efficiency: ' +str(round(maxEta, 4)))
print('Minimum efficiency: ' +str(round(minEta,4)))
print('Valid cases achieved: ' +str(countTrue) + '/'  +str(totalCases) )

# ------------- PLOTTING ------------ 

x = (ZBarr)                     # SAME FOR ALL COUNTOURS
y = np.rad2deg((beta2bArr))     # SAME FOR ALL COUNTOURS
X, Y = np.meshgrid(x, y)  
lvls = 20
colorTheme = 'Reds'

fig, axs21 = plt.subplots(3, 3)
fig.set_figwidth(15)
fig.set_figheight(15)
fig.tight_layout(pad=7.0)
fig.suptitle(r'Flow properties for proposed  $\eta$ = ' + str(etaStage0) +  r' ,  $\lambda _{2}$ = ' + str(lambda2) + ' and Tolerance set to ' +str(iterTol*100) + '% .', y=0.98, x=0.38)
fig.subplots_adjust(top=0.9, bottom=0.09)

# ---------------- Velocity plot ---------------
i=1
j=1
Z = U2Mat
con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
fig.colorbar(con, ax=axs21[i, j])
axs21[i, j].invert_yaxis()
axs21[i, j].set_xticks(x[1::4])
axs21[i, j].set_xticklabels(x[1::4], fontsize=10)
axs21[i, j].set_yticks(np.arange(betamax, betamin , 4))
axs21[i, j].set_yticklabels(np.arange(betamax, betamin , 4), fontsize=10)
axs21[i, j].set_xlabel(r'Blade number $Z_B$ ' , fontsize=12)
axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
axs21[i, j].set_title('Blade exit velocity [m/s]' , fontsize=12)
axs21[i, j].grid()

# ---------------- Flow outlet Velocity plot ---------------
i=1
j=2
Z = c2Mat
con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
fig.colorbar(con, ax=axs21[i, j])
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
fig.colorbar(con, ax=axs21[i, j])
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
fig.colorbar(con, ax=axs21[i, j])
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
fig.colorbar(con, ax=axs21[i, j])
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
cb = fig.colorbar(con, ax=axs21[i, j])
axs21[i, j].invert_yaxis()
axs21[i, j].set_xticks(x[1::4])
axs21[i, j].set_xticklabels(x[1::4], fontsize=10)
axs21[i, j].set_yticks(np.arange(betamax, betamin, 4)) #, fontsize=10)
axs21[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
axs21[i, j].set_title('Exit blade Mach number [-] ' , fontsize=12)
axs21[i, j].grid()


fig, axs2 = plt.subplots(3, 3)
fig.set_figwidth(15)
fig.set_figheight(15)
fig.suptitle(r'Compressor design parameters for proposed $\eta$ = ' + str(etaStage0) +  r' , $\lambda _{2}$ = ' + str(lambda2) + ' and Tolerance set to ' +str(iterTol*100) + '% .', y=0.98, x=0.4)
fig.tight_layout(pad=7.0)
fig.subplots_adjust(top=0.9, bottom=0.09)


# -------------- RPM plot -------------
i=2
j=2
Z = NcritArr                  
con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
fig.colorbar(con, ax=axs2[i, j])
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
cb = fig.colorbar(con, ax=axs2[i, j])
axs2[i, j].invert_yaxis()
axs2[i, j].set_xticks(x[1::4])
axs2[i, j].set_xticklabels(x[1::4], fontsize=10)
axs2[i, j].set_yticks(np.arange(betamax, betamin, 4)) #, fontsize=10)
axs2[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
axs2[i, j].set_ylabel(r' $ \beta _{2B}$ [deg]', fontsize=12)
axs2[i, j].set_title('Inlet hub radius [m]' , fontsize=12)
axs2[i, j].grid()

# -------------- Otlet exit height plot -------------
i=1
j=1
Z = h2Mat                  
con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
cb = fig.colorbar(con, ax=axs2[i, j])
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
cb = fig.colorbar(con, ax=axs2[i, j])
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
fig.colorbar(con, ax=axs2[i, j])
axs2[i, j].invert_yaxis()
axs2[i, j].set_xticks(x[1::4])
axs2[i, j].set_xticklabels(x[1::4], fontsize=10)
axs2[i, j].set_yticks(np.arange(betamax, betamin, 4))
axs2[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
axs2[i, j].set_xlabel(r'Blade number $Z_B$ [deg]', fontsize=12 )
axs2[i, j].set_ylabel(r' $ \beta _{2B}$', fontsize=12)
axs2[i, j].set_title(r'Impeller exit radius [m]', fontsize=12)
axs2[i, j].grid()



print(countcrit)


plt.show(block=True)

