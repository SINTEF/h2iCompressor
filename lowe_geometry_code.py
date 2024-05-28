"""
The following Python code takes a set of variables as inputs and calculates the geometry of an inducer, 
impeller and diffuser for a centrifugal compressor.


Reference: Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244


Author: Martin Spillum Grønli (SINTEF Energy Research, 2024)
"""


### Import-------------------------------------------------------------------------------------------
import math
import numpy as np
import matplotlib.pyplot as plt


### Variables-------------------------------------------------------------------------------------------
P00 = 101300        # Inlet stagnation pressure [Pa]
T00 = 293           # Inlet stagnation temperature [K]
mdot = 0.4          # Mass flow rate [kg/s]
Ndes = 120000       # Rotational speed [rpm]
alpha1 = 0          # Absolute inlet velocity angle [degrees]
rh1 = 0.01          # Hub radius [mm]
B1 = 0.04           # Boundary layer blockage [-]
Cp = 1005           # Specific heat at constant pressure [kJ/kg/K]
k = 1.4             # Ratio of specific heats [-]
R = 287             # Air gas constant [J/kg/K]
Cm1i = np.arange(100, 300.5, 0.5)       # Absolute meridional velocity [m/s]
Pr = 4.8            # Pressure ratio [-]
lambda2 = 2         # Exit swirl parameter
beta2b = - 40       # Exit relative direction [degrees]
eta = 0.6           # Stage efficiency [-]
ZB = 16             # Number of Blades [-]
AR = 2.5            # Area ratio of the diffuser [-]
T00i = np.full(len(Cm1i), 293)      # Inlet stagnation temperature [K]
D2 = 0.105          # Chosen exit diameter [m]


### Inducer calculation----------------------------------------------------------------------------------
"""
For each value of Cm1 in Cm1i, we calculate the corresponding values of C1, T1, M1, P1, rho1, A1, rt1, U1t, W1t, Cm1 and beta1.
The value of Cm1 will affect the values of U1t and Ctheta1 (through ex. rt1), which in turn affects the value of W1t.
By plotting W1t against Cm1, we then find the value of Cm1 that minimises W1t.
"""
Ctheta1i = Cm1i * math.tan(math.radians(alpha1))    # Inlet absolute tangential velocity [degrees]
C1i = (Ctheta1i ** 2 + Cm1i ** 2) ** 0.5            # Absolute velocity [m/s]
T1i = T00i - (C1i ** 2) / (2 * Cp)                  # Inlet temperature [K]
M1i = C1i / ((k * R * T1i) ** 0.5)                  # Inlet Mach number [-]
P1i = P00 * (T1i / T00) ** (k / (k - 1))            # Inlet pressure [Pa]
rho1i = P1i / (R * T1i)                             # Inlet density of air [kg/m^3]
A1i = mdot / (rho1i * Cm1i * (1 - B1))              # Inlet flow area [m^2]
#rt1i = (A1i / (math.pi * (1 - (hubtip) ** 2))) ** 0.5      # Inlet tip radius [m]
rt1i = (A1i / math.pi + rh1 ** 2) ** 0.5            # Inlet tip radius [m]
U1ti = 2 * math.pi * rt1i * Ndes / 60               # Inlet blade tip speed [m/s]
W1ti = (Cm1i ** 2 + ((U1ti ** 2) - (Ctheta1i ** 2))) ** 0.5     # Inlet relative velocity [m/s]     # MSG: Should be (Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5?
#W1ti = (Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5
Ctheta1 = Ctheta1i[np.argmin(W1ti)]     # Get the value of Ctheta1 for minimal value of W1t
C1 = C1i[np.argmin(W1ti)]
T1 = T1i[np.argmin(W1ti)]
M1 = M1i[np.argmin(W1ti)]
P1 = P1i[np.argmin(W1ti)]
rho1 = rho1i[np.argmin(W1ti)]
A1 = A1i[np.argmin(W1ti)]
rt1 = rt1i[np.argmin(W1ti)]
U1t = U1ti[np.argmin(W1ti)]
W1t = W1ti[np.argmin(W1ti)]
Cm1 = Cm1i[np.argmin(W1ti)]
beta1 = math.degrees(math.atan((U1t - Ctheta1) / Cm1))      # Inlet relative velocity [degrees]

def W1tminimise(Cm1i, U1ti, Ctheta1i):      # MSG: Rename function (not a minimization function)
    return (Cm1i ** 2 + ((U1ti ** 2) - (Ctheta1i ** 2))) ** 0.5     # MSG: Should be (Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5?
    #return (Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5

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

W1tmin = W1tminimise(Cm1i, U1ti, Ctheta1i)      # Get W1t as a function of Cm1
plt.figure()
plt.plot(Cm1i, W1tmin)
plt.xlabel("Cm1 (m/s)")
plt.ylabel("W1t (m/s)")
plt.title("Minimisation of W1t")
plt.show()


### Impeller calculation---------------------------------------------------------------------------------
sigma = 1 - (math.sqrt(math.cos(math.radians(beta2b))) / (ZB ** 0.7))       # Slip factor [-]   MSG: Can formula be adjusted? This is probably a correlation, see e.g. Eqs. (59)-(60) in HYDROGENi memo D2.3.5 Current understanding and knowledge gaps for operation of hydrogen gas export compression systems.
mu = sigma * lambda2 / (lambda2 - math.tan(math.radians(beta2b)))           # Work input coefficient [-]    MSG: Can formula be adjusted? 
hx = ((k * R * T00) / (k - 1)) * (Pr ** ((k - 1) / k) - 1)                  # Specific enthalpy [kJ/kg/K]
Wx = hx / eta                               # Specific work [kJ/kg/K]
T02m = T00 + Wx * (k - 1) / (k * R)         # Stagnation exit temperature [K]
P02m = P00 * (((k - 1) * Wx * eta / (k * R * T00)) + 1) ** (k / (k - 1))    # Exit Stagnation Pressure [Pa]
#U2 = ((U1t * Ctheta1 + Wx) / mu) ** 0.5    # Exit Blade Speed [m/s]
#D2 = 60 * U2 / (math.pi * Ndes)            # Exit Diameter [m]
U2 = D2 * math.pi * Ndes / 60               # Exit blade speed [m]
Ctheta2m = mu * U2                          # Absolute tangential exit velocity [m/s]
Cm2m = Ctheta2m / lambda2                   # Absolute meridional exit velocity [m/s]
T2m = T02m - (k - 1) / (2 * k * R) * (Ctheta2m ** 2 + Cm2m ** 2)            # Exit temperature [K]
P2m = P02m / ((T02m / T2m) ** (k / (k - 1)))                                # Exit pressure [Pa]
rho2m = P2m / (T2m * R)                     # Exit density [kg/m^3]
A2 = mdot / (rho2m * Cm2m)                  # Exit area [m^2]
b2 = A2 / (math.pi * D2)                    # Depth of impeller exit [m]
C2 = (Ctheta2m ** 2 + Cm2m ** 2) ** 0.5     # Absolute exit velocity [m/s]


### Diffuser Calculation----------------------------------------------------------------------------------
etad = 0.85                                 # Estimated diffuser efficiency
CpDi = 1 - (AR ** - 2)                      # Ideal pressure recovery coefficient
CpD = etad * CpDi                           # Pressure recovery coefficient
P3 = P2m + CpD * (P02m - P2m)               # Diffuser exit static pressure [Pa]
C3 = C2 / AR                                # Diffuser exit absolute velocity [m/s]
P03 = P3 + 0.5 * rho2m * C3 ** 2            # Diffuser exit stagnation pressure [Pa]


### Overall performance----------------------------------------------------------------------------------
etaiterate = ((P03 / P00) ** ((k - 1) / k) - 1) / ((T02m / T00) - 1)        # Iterative stage efficiency [-]
Prest = ((etaiterate * U2 ** 2 * mu) / (Cp * T1) + 1) ** (k / (k - 1))      # Estimate of the pressure ratio


### Code Outputs------------------------------------------------------------------------------------------
print("Pressure ratio =", P03 / P00)
print("Pressure ratio error =", Pr - P03 / P00)
print("Efficiency =", etaiterate)
print("Efficiency error =", abs(etaiterate - eta))
print("Work =", Wx * mdot, "Watts")
print("Mass flow rate =", mdot, "kg/s")
print("Exducer diameter =", D2 * 1000, "mm")
print("Inducer Diameter =", rt1 * 2000, "mm")
print("Exit blade speed =", U2)
if rt1 / (0.5 * D2) < math.exp(- 8.16 * math.cos(beta2b) / ZB):     # Check for sigma validity
    print("sigma valid")
else:
    print("sigma invalid!")