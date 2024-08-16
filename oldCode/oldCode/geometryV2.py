"""
The following Python code takes a set of variables as inputs and calculates the geometry of an inducer, 
impeller and diffuser for a centrifugal compressor.


Reference: Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244


Author: Petter Resell (SINTEF Energy Research, 2024)
"""


### Import----------------------------g---------------------------------------------------------------
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#plt.close('all')

# Set plot parameters------------------------------------------------------------------------------
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams.update({'font.size': 15})


### Variables-------------------------------------------------------------------------------------------


# Inlet flow parameters
mdot = 20                         # Mass flow rate [kg/s]   
# mdot = ( (300* 10**3 )/( 60**2 ) )/6    #  => todo = 

Ndes = 120000                           # Rotational speed [rpm]
Ndes = 60000
Ndes = 10000                          # Rotational speed [rpm]
# fra memo figur 

Cm1i = np.arange(10, 600.5, 0.5)       # Inlet Absolute meridional velocity [m/s]

# Fluid properties  
CpAir = 1006                            # Cp air [J/kg/K], 
MolarMassAir = 0.02897                  # Molecular weight of air [kg/mol]. 
CpH2 = 14310                            # Engineering toolbox Cp hydrogen, [J/kg/K]                                       
MolarMassH2 = 2.01568 * 10**-3          # Molecular weight of H2 [kg/mol]. 
Cp = CpH2                          # Engineering toolbox Cp hydrogen, [J/kg/K]                                       
MolarMass = MolarMassH2          # Molecular weight of H2 [kg/mol]. 

k = 1.41
                               # Ratio of specific heats, isentropic exponent [-],    apply to both O2 and H2                             
R_uni = 8.314                         # Universal gas constant [J/mol/K]
R = R_uni / MolarMass                 # Specific gas constant [J/kg /K]

# Impeller exit parameters
lambda2 = 2                 # Exit swirl parameter                                          , check out
lambda20 = lambda2          # Exit swirl parameter, used for iteration
beta2b = -40                # Exit relative direction [deg], forward leaning blade -40deg. Typically on interval [-10, -45]                    25
beta2b0 = beta2b            # Exit relative direction [deg], used for iteration                             
etaStage = 0.6              # Isentropic Stage efficiency [-]                                          
etaiterate = etaStage       # Isentropic stage efficiency used for iteration [-]

# Generall compressor parameters
# inlet conditions
P00 = 30 * (10**2) * (10**3)                    # Inlet stagnation pressure [Pa]
T00 = 293                                       # Inlet stagnation temperature [K]
Pr = 1.2                   # Pressure ratio [-]
# Pr = 2

alpha1 = 0                                      # Absolute inlet velocity angle [degrees]
rh1 = 0.1                                       # Hub radius [m]
B1 = 0.04                                       # Boundary layer blockage [-],        substance viscocity dependant
ZB = 5                                          # Number of Blades, 8 splitter blades & 8 full blades. Must be even number [-]        , typisk mellom 10 og 20 blades pressure ratios 2.5->6
ZB0 = ZB                                        # Number of blades, used for iteration
AR = 2.5                                        # inlet/outlet area ratio of the diffuser [-]         
T00i = np.full(len(Cm1i), 293)                  # Inlet stagnation temperature [K]
D2 = 0.105                                      # Chosen impeller exit diameter [m]         
# D2 = 0.43                                     # Chosen impeller exit diameter [m]         

omega = 2 * math.pi * Ndes/60       # Rotational velocity of blades [rad/sec]


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
rt1i = (A1i / math.pi + rh1 ** 2) ** 0.5                    # Inlet tip radius [m]                           eqn (52) rapport
U1ti = 2 * math.pi * rt1i * Ndes / 60                       # Inlet blade tip speed [m/s]                    eqn (60) rapport

#Todo
""" which of the two equations underneith is used doesnt matter when inlet flow angle is alpha1 zero. Ctheta1 is then zero-array """
#W1ti = (Cm1i ** 2 + ((U1ti ** 2) - (Ctheta1i ** 2))) ** 0.5     # Inlet relative velocity [m/s]     # Todo MSG: Should be (Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5?
W1ti = (Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5        # Todo    : enig i denne linjen, ikke den over

### kvalitetssikret hittil => fortsett onsdag

# Todo: Hvorfor ønsker man å minimere W1ti??? => efficiency

"""Problemstilling: rt1 > rt2, må iterere å finne den indexen som gir minst relativ hastighet samtidig som r2t>rt1
    EDIT: rt1 er alltid større enn rt2  => PROBLEM"""

rt1 = D2*0.5                    # initializing value, is changed immediately
minimizeCount = 1               # Counting iterations
W1tiIterate = np.copy(W1ti)     # W1tIterate = W1ti made w1ti change inside loop also

# todo: make this a function
# Not implemented at this point in time
while rt1 >= 0.5*D2 and minimizeCount <= len(W1tiIterate):
    indice_min = np.argmin(W1tiIterate)             # np.argmin returns indice of minimum value
    rt1 = rt1i[indice_min]
    if rt1 >= D2*0.5:
        W1tiIterate[indice_min] = max(W1tiIterate)
  
    if minimizeCount == len(W1tiIterate):
        indice_min = np.argmin(W1ti)
        rt1 = rt1i[indice_min]
    minimizeCount +=1

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

""" Unnessecary to display:
print("W1t is minimal for:")
print("\tW1t =", W1t)
print("\tC1 =", C1)
print("\tCm1 =", Cm1)
print("\tCtheta1 =", Ctheta1)
print("\tU1t =", U1t)
print("\tbeta1 =", beta1)
print("\tT1 =", T1)
print("\tM1 =", M1)
print("\tP1 =", P1)
print("\trho1 =", rho1)
print("\tA1 =", A1)
print("\trt1 =", rt1)
"""

W1tIter = W1tIterate(Cm1i, U1ti, Ctheta1i)      # Get W1t as a function of Cm1
plt.figure()
plt.plot(Cm1i, W1tIter)
plt.xlabel("Cm1 (m/s)")
plt.ylabel("W1t (m/s)")
plt.grid()
plt.title("Minimisation of W1t")

# -------------------------------------------------------------------------------------------------------------

# Iterate-function decides which  variable to change when iterating. Change this by changing default arguments a1, a2 & a3 
def iterateInnerLoop(beta2bArg = beta2b, ZbArg = ZB, lambda2Arg = lambda2, a1 = False, a2 = True, a3 = False, countInnerLoopCheck=1):
        if a1 == True and countInnerLoopCheck>1:
            beta2bArg += 0.5
        if a2 == True  and countInnerLoopCheck>1:
            ZbArg += 1
        if a3 == True and countInnerLoopCheck>1:
            lambda2Arg += 0        # Todo: find fitting adjustment

        # Checking edge case

        if ZbArg <= 0 or ZbArg > 100:
            print('\n')
            raise Exception("Cant have zero, negatve or unpractically large number of impeller blades: " + str(ZbArg))
        
        if beta2b >= 0 or beta2b < -80:
            raise Exception('beta2b unvalid: ' + str(beta2bArg))
                  
        return [beta2bArg, ZbArg, lambda2Arg]

### Impeller calculation---------------------------------------------------------------------------------
"""The isentropic enthalpy demand is used to find the approximate required work. Furthermore the  Work 
is used to find the outlet stagnation temperature which in turn is used to find the outlet stagnation pressure 
from isentropic relations. Further, geometrical relations are applied. - Petter """

# Defining iteration variables to keep track of looping
countOuterLoop = 1                  # outer loop iteration count
countInnerLoop = 1                  # Inner loop iteration count
checkOuterLoop = False              # T/F outer loop check
PressureTestOuterLoop = 1           # Testvariable outer loop


"""If the pressure estimate Prest initially over-estimates the pressure ratio by more than 5%, then the the error will only grow for increasing efficiency
        since the pressure ratio increase for increasing efficiency. Hence its necessary to check if the pressure estimate is an 
            under- or over-estimate before adjusting the efficiency eta. From then the efficiency can be either increased or decreased. """

def pressureOverUnderEstimate(count, outerTest):
    if count == 1 or (countOuterLoop > 1 and outerTest > 0):
        etaST = etaiterate                  # todo: need better way of iterating etaStage
    elif count > 1 and outerTest < 0:
        etaST = etaStage + 0.01
    return etaST



""" ------------------------------------ Beginning iteration ------------------------------------ """
print('\n')

iterTolOuterLoop = 0.025            # Outer loop condition tolerance

dh0s = ( (k * R * T00) / (k - 1) ) * (Pr **((k - 1)/ k) - 1)                     # Specific isentropic compression enthalpy, constant value [kJ/kg/K]

while checkOuterLoop == False or countOuterLoop == 1:

    checkInnerLoop = False
    beta2b, ZB, lambda2 = iterateInnerLoop( beta2bArg=beta2b0, ZbArg=ZB0, lambda2Arg=lambda20, a1=0, a2=0, a3=0 )

    print("Outer loop iteration number ", countOuterLoop)

    while checkInnerLoop == False:
        innerloopRHS = math.exp(-8.16*math.cos(np.deg2rad(beta2b))/ZB)          # Inner loop RHS, exponential expression

        # Wiesner slip factor
        """ Slip factor [-]   MSG: Can formula be adjusted? This is probably a correlation, see e.g. Eqs. (59)-(60) in HYDROGENi memo
        D2.3.5 Current understanding and knowledge gaps for operation of hydrogen gas export compression systems.
        As per now only a function of exit angle and number of blades ZB
        Maybe try to implement:
        - Stodolas equation
        - Stainzs equation
        but both depend on work coefficient found later
        - Balje's formula

        => https://en.wikipedia.org/wiki/Slip_factor
        """
        Wx = dh0s / etaStage                                                             # Specific work [kJ/kg/K]
        sigma = 1 - (math.sqrt(math.cos(math.radians(beta2b))) / (ZB ** 0.7))            # slip factor
        T02m = T00 + Wx * (k - 1) / (k * R)                                              # Stagnation exit temperature [K]     , from dh=cp*Dt   (k - 1) / (k * R) = 1/cp
        P02m = P00 * ((etaStage * Wx * (k - 1) / (k * R * T00)) + 1) ** (k / (k - 1))    # Exit Stagnation Pressure [Pa]       , from (... todo ...)

        mu = sigma * lambda2 / (lambda2 - math.tan(math.radians(beta2b)))                # Work input coefficient [-]    MSG: Can formula be adjusted? tan(-x=-tan(x)) 
        # Todo: Mener etaStage i uttrykket over skal bort burde være:
        # P02m = P00*(T02m/T00)**(k/(k-1))

        U2 = ((U1t * Ctheta1 + Wx) / mu) ** 0.5    # Exit Blade Speed [m/s]
        D2 = 60 * U2 / (math.pi * Ndes)              # Exit Diameter [m]
        # U2 = D2 * math.pi * Ndes / 60
        # !Dont touch!
        Ctheta2m = mu * U2                                # Absolute tangential exit velocity [m/s]         from work coefficient
        Cm2m = Ctheta2m / lambda2                         # Absolute meridional exit velocity [m/s]         
        # Cm2m = (U2 - Ctheta2m)/math.tan(beta2b)         # Absolute meridional exit velocity [m/s] alternative 
        phi2 = Cm2m/U2                                    # Flow coeff. [-]
        # sigma = 1 - (math.sqrt(math.cos(math.radians(beta2b))) / ( (ZB ** 0.7)*(1 - phi2 * math.tan(beta2b)) ))  # Alternative for slip factor iterere til kopnvergens????

        C2 = (Ctheta2m ** 2 + Cm2m ** 2) ** 0.5           # Absolute exit velocity [m/s]                    from pythagoras 
        T2m = T02m - (k - 1) / (2 * k * R) * (C2 **2)     # Exit temperature [K]                            from stagnation temperature
        P2m = P02m / ((T02m / T2m) ** (k / (k - 1)))      # Exit pressure [Pa]                              from stagnation pressure
        rho2m = P2m / (T2m * R)                           # Exit density [kg/m^3]                           from ideal gas law
        A2 = mdot / (rho2m * Cm2m)                        # Exit area [m^2]                                 from continuity
        b2 = A2 / (math.pi * D2)                          # Depth of impeller exit [m]                      from geometry

        
        # Displaying results
        print(' Etaiterate: ', etaStage , ' Beta2: ', beta2b, ' Zb: ', ZB, '              Slip factor: ', sigma, ' r1/r2: ', rt1/(0.5*D2), '    RHS inner condition: ', innerloopRHS)
        countInnerLoop +=1

        r1tDivr2t = rt1/(0.5*D2)   
        # Checking inner loop condition
        if r1tDivr2t >= innerloopRHS and (ZB <= 30 or r1tDivr2t) :
            checkInnerLoop = False
            beta2b, ZB, lambda2 = iterateInnerLoop(beta2b, ZB, lambda2, countInnerLoopCheck=countInnerLoop)     # Updating vaqriables for next iteration
        elif r1tDivr2t > 0.7 or ZB > 30:
             checkInnerLoop = False
             break
        else:
            checkInnerLoop = True

                                         # LHS of inner loop condition, ratio of radii
        # Checking edge case
        if beta2b >= 0 or beta2b < -80:
            raise Exception('beta2b greater than zero')




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
    if abs(PressureTestOuterLoop) > iterTolOuterLoop or checkInnerLoop == False:
        checkOuterLoop = False
    else:
        checkOuterLoop = True


    """If the pressure estimate Prest initially over-estimates the pressure ratio by more than 5%, then the the error will only grow for increasing efficiency
        since the pressure ratio increase for increasing efficiency. 6
        Hence its neeeded to check if the pressure estimate is an under- or over-estimate before adjusting the efficiency eta. """


    etaStage = pressureOverUnderEstimate(countOuterLoop, PressureTestOuterLoop)
    if etaStage >= 1:
        raise Exception("efficiency cant be greater or equal to 1 => CHECK TOLERANCE") 
    countOuterLoop +=1

    print(' Test outer loop: ', round(abs(PressureTestOuterLoop), 10), '                    New stage eta: ', etaStage)
    print('\n')


### Code Outputs------------------------------------------------------------------------------------------
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
print("Inducer Diameter =", round(rt1 * 2000, digits), "mm")
print("Exit blade speed =", round(U2, digits))
print("Number of iterations of outer loop: ", int(countOuterLoop-1))
"""

# if rt1 / (0.5 * D2) < math.exp(- 8.16 * math.cos(beta2b) / ZB):     # Check for sigma validity
#     print("sigma valid")
# else:
#     print("sigma invalid!")

data1 = [round(P03 / P00, digits), round(Pr - P03 / P00, digits), round(etaiterate, digits), round( (abs(etaiterate - etaStage)), digits), 
          round(Wx * mdot, digits-3), round(mdot, digits), round(D2 * 1000, digits),round(rt1 * 2000, digits),int(round(U2, digits)), P00* 10**-5, P03* 10**-5 ]
col1 = ["Pressure ratio", "Pressure ratio error", "Efficiency", "Efficiency error", "Work", "Mass flow rate","Exducer diameter", 
        "Inducer Diameter", "Exit blade speed", "Inlet Pressure", "Outlet Pressure"]
units = ['-', '-', '-', '-', 'Watt', 'kg/s', 'mm', 'mm', '?', '[bar]', '[bar]']
data = {"*Variable*    ":col1, "      *Value*":data1, "      *Unit*":units}
df = pd.DataFrame(data)
table = df.to_string(index=False)
print(table)
print('\n')

plt.show()