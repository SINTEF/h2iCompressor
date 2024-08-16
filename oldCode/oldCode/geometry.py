"""
The following Python code takes a set of variables as inputs and calculates the geometry of an inducer, 
impeller and diffuser for a centrifugal compressor.


Reference: Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244


Author: Petter Resell (SINTEF Energy Research, 2024)
"""
# Todo:
# - utdyp forklaring av koden, og kommenter hva som skjer
# - 
# - 
# - 
# - 

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

# Todo: parametriser ulike caser i excel og kjør loop

P00 = 101300        # Inlet stagnation pressure [Pa]
T00 = 293           # Inlet stagnation temperature [K]
mdot = 0.4          # Mass flow rate [kg/s], 300tonn/time  , /4     todo 300tonn/time = 300*10^3 kg/time = 300* 10**3 / 60**2
mdot = ( (300* 10**3 )/( 60**2 ) )/6
Ndes = 120000       # Rotational speed [rpm]
alpha1 = 0          # Absolute inlet velocity angle [degrees]
rh1 = 0.01          # Hub radius [m]
B1 = 0.04           # Boundary layer blockage [-]                                                           substance viscocity dependant
Cp = 1006           # Specific heat at constant pressure [kJ/kg/K],                                         substance dependant
k = 1.4             # Ratio of specific heats, isentropic exponent [-],                                     substance dependant
#MolarMass = 0.02897         # Molecular weight of air [kg/mol]. 
MolarMass = 2.01568* 10**-3          # Molecular weight of H2 [kg/mol]. 
R_uni = 8.314       # Universal gas constant [J/mol/K]
R = R_uni / MolarMass       # Specific gas constant [J/kg/K]
Cm1i = np.arange(100, 300.5, 0.5)       # Inlet Absolute meridional velocity [m/s]
Pr = 4.8            # Pressure ratio [-]
lambda2 = 2         # Exit swirl parameter                                          , sjekk ut
lambda20 = lambda2
beta2b0 = -40   # Exit relative direction [degrees]                             , forward leaning blades (?) -40 base
beta2b = beta2b0       # Exit relative direction [degrees], blade exit angle                             , forward leaning blades (?)
etaStage = 0.6           # Isentropic Stage efficiency [-]                                          , polytropic?
etaiterate = etaStage
etaMax = etaStage
ZB = 16             # Number of Blades, 8 splitter blades, 8 full blades [-]        , have to keep number even to represent splitter blades + full blades 
ZB0 = ZB
AR = 2.5            # Area ratio of the diffuser [-]                                , inlet/outlet area
T00i = np.full(len(Cm1i), 293)      # Inlet stagnation temperature [K]
D2 = 0.105          # Chosen exit diameter [m]                                      , Impeller exit???




### Inducer calculation----------------------------------------------------------------------------------
"""
For each value of Cm1 in Cm1i, we calculate the corresponding values of C1, T1, M1, P1, rho1, A1, rt1, U1t, W1t, Cm1(?) and beta1.
The value of Cm1 will affect the values of U1t and Ctheta1 (through ex. rt1), which in turn affects the value of W1t.
By plotting W1t against Cm1, we then find the value of Cm1 that minimises W1t.

rt1 inlet tip radius ???

"""
                                                                # what                                       # from
Ctheta1i = Cm1i * math.tan(math.radians(alpha1))            # Inlet absolute tangential velocity [degrees]
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
W1ti = (Cm1i ** 2 + ((U1ti ** 2) - (Ctheta1i ** 2))) ** 0.5     # Inlet relative velocity [m/s]     # Todo MSG: Should be (Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5?
#W1ti = (Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5        # Todo    : enig i denne linjen, ikke den over


# Todo: Hvorfor ønsker man å minimere W1ti???
indice_min = np.argmin(W1ti)            # np.argmin returns indice of minimum value
# print('max: ', max(W1ti))
# print('min: ', min(W1ti))

Ctheta1 = Ctheta1i[indice_min]     # Get the value of Ctheta1 for minimal value of W1t
C1 = C1i[indice_min]        
T1 = T1i[indice_min]
M1 = M1i[indice_min]
P1 = P1i[indice_min]
rho1 = rho1i[indice_min]
A1 = A1i[indice_min]
rt1 = rt1i[indice_min]
U1t = U1ti[indice_min]
W1t = W1ti[indice_min]
Cm1 = Cm1i[indice_min]
beta1 = math.degrees(math.atan((U1t - Ctheta1) / Cm1))      # Inlet relative velocity [degrees]

def W1tIterate(Cm1i, U1ti, Ctheta1i):                              # Todo MSG: Rename function (not a minimization function)
    #return (Cm1i ** 2 + ((U1ti ** 2) - (Ctheta1i ** 2))) ** 0.5     # Todo MSG: Should be (Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5?
    return (Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5

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

W1tIter = W1tIterate(Cm1i, U1ti, Ctheta1i)      # Get W1t as a function of Cm1
plt.figure()
plt.plot(Cm1i, W1tIter)
plt.xlabel("Cm1 (m/s)")
plt.ylabel("W1t (m/s)")
plt.grid()
plt.title("Minimisation of W1t")

""""Iterate-function decides which """
globVar=1
def iterateInnerLoop(beta2bArg = beta2b, ZbArg = ZB, lambda2Arg = lambda2, a1 = 0, a2 = 1, a3 = 0, countInnerLoopCheck=1):
        if a1 == 1 and countInnerLoopCheck>1:
            beta2bArg -= 0.5
        if a2 == 1 and countInnerLoopCheck>1:
            ZbArg -= 2
        if a3 == 1 and countInnerLoopCheck>1:
            lambda2Arg += 0        # Todo: find fitting adjustment
        if ZbArg <= 0:
            raise Exception("Cant have zero or negatve number of impeller blades")
        return [beta2bArg, ZbArg, lambda2Arg]

### Impeller calculation---------------------------------------------------------------------------------
"""The isentropic enthalpy demand is used to find the approximate required work. Furthermore the  Work 
is used to find the outlet stagnation temperature which in turn is used to find the outlet stagnation pressure 
from isentropic relations. Further, geometrical relations are applied. - Petter """


countOuterLoop = 1
countInnerLoop = 1
iterTolOuterLoop = 0.05
checkInnerLoop = False
checkOuterLoop = False
testOuterLoop = 1 


"""If the pressure estimate Prest initially over-estimates the pressure ratio by more than 5%, then the the error will only grow for increasing efficiency
        since the pressure ratio increase for increasing efficiency. Hence its necessary to check if the pressure estimate is an 
            under- or over-estimate before adjusting the efficiency eta. From then the efficiency can be either increased or decreased. """

def pressureOverUnderEstimate(count, outerTest):
    etaST = 0
    if count == 1 or (countOuterLoop > 1 and outerTest > 0):
        etaST = etaiterate        # todo: need better way of iterating etaStage
    elif count > 1 and outerTest < 0:
        etaST = etaiterate + 0.005
    return etaST

    



""" ------------------------------------ Beginning iteration ------------------------------------ """

while checkOuterLoop == False:
    # Todo implement inner loop
    # initialize(check2init=checkOuterLoop)
    checkInnerLoop = False
    beta2b, ZB, lambda2 = iterateInnerLoop( beta2bArg=beta2b0, ZbArg=ZB0, lambda2Arg=lambda20, a1=0, a2=0, a3=0 )

    print("Outer loop iteration number ", countOuterLoop)

    while checkInnerLoop == False:
        
        # Wiesner slip factor
        sigma = 1 - (math.sqrt(math.cos(math.radians(beta2b))) / (ZB ** 0.7))       # slip factor
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

        mu = sigma * lambda2 / (lambda2 - math.tan(math.radians(beta2b)))            # Work input coefficient [-]    MSG: Can formula be adjusted? 
        dhx = ((k * R * T00) / (k - 1)) * (Pr ** ((k - 1) / k) - 1)                  # Specific isentropic compression enthalpy [kJ/kg/K]
        Wx = dhx / etaStage                               # Specific work [kJ/kg/K]
        T02m = T00 + Wx * (k - 1) / (k * R)         # Stagnation exit temperature [K]                                           from dh=cp*Dt   (k - 1) / (k * R) = 1/cp
        P02m = P00 * ((etaStage * Wx * (k - 1) / (k * R * T00)) + 1) ** (k / (k - 1))    # Exit Stagnation Pressure [Pa]             from (... todo ...)
        # Todo: Mener etaStage i uttrykket over skal bort burde være:
        # P02m = P00*(T02m/T00)**(k/(k-1))

        #U2 = ((U1t * Ctheta1 + Wx) / mu) ** 0.5    # Exit Blade Speed [m/s]
        #D2 = 60 * U2 / (math.pi * Ndes)            # Exit Diameter [m]

        U2 = D2 * math.pi * Ndes / 60               # Exit blade speed [m]                            from rotational velocity relations
        Ctheta2m = mu * U2                          # Absolute tangential exit velocity [m/s]         from work coefficient
        Cm2m = Ctheta2m / lambda2                   # Absolute meridional exit velocity [m/s]         from geometry of velocity triangle
        C2 = (Ctheta2m ** 2 + Cm2m ** 2) ** 0.5     # Absolute exit velocity [m/s]                    from pythagoras 
        T2m = T02m - (k - 1) / (2 * k * R) * (C2 **2)     # Exit temperature [K]                      from stagnation temperature
        P2m = P02m / ((T02m / T2m) ** (k / (k - 1)))      # Exit pressure [Pa]                        from stagnation pressure
        rho2m = P2m / (T2m * R)                     # Exit density [kg/m^3]                           from ideal gas law
        A2 = mdot / (rho2m * Cm2m)                  # Exit area [m^2]                                 from continuity
        b2 = A2 / (math.pi * D2)                    # Depth of impeller exit [m]                      from geometry

        # Check1-variable to be greater than 1;
        testInnerLoop = ( 0.5*D2/rt1)*math.exp(-8.16*math.cos(beta2b)/ZB)
        if testInnerLoop <= 1:
            checkInnerLoop = False
            beta2b, ZB, lambda2 = iterateInnerLoop(beta2b, ZB, lambda2, countInnerLoopCheck=countInnerLoop)
        else:
            checkInnerLoop = True

        print('Test inner loop: ', round(testInnerLoop, 5), ' Etaiterate: ', etaStage , ' Beta2: ', beta2b, ' Zb: ', ZB, '              Slip factor: ', sigma, ' r1/r2: ', rt1/(0.5*D2), 'RHS inner condition: ', (rt1/(0.5*D2)) * testInnerLoop)
        
        countInnerLoop +=1
        
        
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


       
    testOuterLoop = (Prest-Pr)/Pr
    if abs(testOuterLoop) > iterTolOuterLoop:
        checkOuterLoop = False
    else:
        checkOuterLoop = True


    """If the pressure estimate Prest initially over-estimates the pressure ratio by more than 5%, then the the error will only grow for increasing efficiency
        since the pressure ratio increase for increasing efficiency. 
        Hence its neeeded to check if the pressure estimate is an under- or over-estimate before adjusting the efficiency eta. """


    etaStage = pressureOverUnderEstimate(countOuterLoop, testOuterLoop)
    if etaStage >= 1:
        raise Exception("efficiency cant be greater or equal to 1") 
    countOuterLoop +=1

    print('Test outer loop: ', round(abs(testOuterLoop), 10), '         New stage eta: ', etaStage)
    print('\n')


### Code Outputs------------------------------------------------------------------------------------------
digits = 5
print('--------------------------- Code Outputs ---------------------------')
print("Pressure ratio =", round(P03 / P00, digits) )
print("Pressure ratio error =", round(Pr - P03 / P00, digits))
print("Efficiency =", round(etaiterate, digits))
print("Efficiency error =", round( (abs(etaiterate - etaStage)), 5) )
print("Work =", round(Wx * mdot, digits-3), "Watts")
print("Mass flow rate =", round(mdot, digits), "kg/s")
print("Exducer diameter =", round(D2 * 1000, digits), "mm")
print("Inducer Diameter =", round(rt1 * 2000, digits), "mm")
print("Exit blade speed =", round(U2, digits))
print("Number of iterations of outer loop: ", int(countOuterLoop-1))
if rt1 / (0.5 * D2) < math.exp(- 8.16 * math.cos(beta2b) / ZB):     # Check for sigma validity
    print("sigma valid")
else:
    print("sigma invalid!")

data1 = [round(P03 / P00, digits), round(Pr - P03 / P00, digits), round(etaiterate, digits), round( (abs(etaiterate - etaStage)), digits), 
          round(Wx * mdot, digits-3), round(mdot, digits), round(D2 * 1000, digits),round(rt1 * 2000, digits),round(U2, digits), countOuterLoop]
col1 = ["Pressure ratio", "Pressure ratio error", "Efficiency", "Efficiency error", "Work", "Mass flow rate","Exducer diameter", 
        "Inducer Diameter", "Exit blade speed","Number of iterations of outer loop",]
units = ['-', '-', '-', '-', 'Watt', 'kg/s', 'mm', 'mm', '?','-' ]
data = {"*Variable*    ":col1, "      *Value*":data1, "      *Unit*":units}
df = pd.DataFrame(data)
table = df.to_string(index=False)
print(table)

plt.show()
