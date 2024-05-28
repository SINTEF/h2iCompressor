"""
The following Python code uses the inducer and impeller geometry from the geometry code
(FinalDesignCode) and predicts the performance of the compressor at off-design points.

Note that the letter B and number inside the brackets in some of the lines of code correspond
to the equations in NASA’s Fortran model (B46, B47…B60 for example). The B denotes
Appendix B in the original paper and the number corresponds to the equation number in
approximate order of use in the code (Galvas, 1973). There is a jump in numbers between B7-
B40 and B72-B84. The missing equations correspond to inlet guide vane losses and vaned
diffuser losses, both of which were omitted for reasons outlined in Section 4.2. Many of the
variables and equations used in this code have not been defined in this thesis. The full set of
equations and a description of the calculation process can be found in Appendix B of NASA’s
Fortran model. A description of each calculation line in the code is provided after the hash
symbol (#). Additionally, a description of each function that has been used is also provided.


Reference: Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244


Author: Martin Spillum Grønli (SINTEF Energy Research, 2024)
"""


### Import---------------------------------------------------------------------------------------------------
from lowe_geometry_code import *
import math
import numpy as np
import matplotlib.pyplot as plt


### Variables------------------------------------------------------------------------------------------------
Ndes = 120e3
No = 1.0            # Off design speed percentage that requires performance prediction [rpm]
tu = 0.002          # Blade thickness [m]
rho0 = 1.1839       # Stagnation density of air [kg/m^3]
curvet1 = 0         # Inducer inlet tip wall curvature [m^-1]
curveh1 = 0         # Inducer inlet hub wall curvature [m^-1]
x = 10              # Streamline angle [degrees] from axial direction
VCR = math.sqrt(2 * k / (k + 1) * R * T00)      # Critical Velocity [m/s]
VOVCR = np.linspace(0.1, 0.68, 30)              # Compressor inlet absolute critical velocity ratio [-]
Cm1h = VOVCR * VCR  # Absolute meridional velocity of the hub [m/s]
kBL = 0.6           # Blading loss coefficient [-]
visc = 1.778e-5     # Air viscosity based on total conditions
kSF = 7.0           # Skin friction coefficient [-]
Cf = 0.01           # Friction coefficient [-]
T00i = np.full(len(VOVCR), 293)


### Calculation of Inlet Velocity Triangles and Compressor Weight Flow (Swirl free)----------------------------
curve1rms = math.sqrt((curvet1 ** 2 + curveh1 ** 2) / 2)        # Root mean square of the inducer inlet hub and tip wall curvature [m^-1]
r1rms = math.sqrt((rt1 ** 2 + rh1 ** 2) / 2)                    # Root mean square of the hub and tip inlet radius [m]

h0 = r1rms - rh1            # (B4) Spacing for numerical integration [-]
h1 = rt1 - r1rms            # (B5) Spacing for numerical integration [-]
Cm1rms = Cm1h * math.exp((h0 / 2) * (curveh1 + curve1rms))      # (B2) Absolute meridional root mean square velocity [m/s]
Cm1t = Cm1h * math.exp(((h0 + h1) / 6) * ((2 - (h1 / h0)) * curveh1 + (((h0 + h1) ** 2) / (h1 * h0)) * curve1rms + (2 - (h0 / h1)) * curvet1))      # (B3) Absolute meridional velocity of the tip [m/s]
Cm1hn = Cm1h * math.cos(math.radians(x))        # (B6) Normal component of the absolute hub velocity [m/s]
Cm1tn = Cm1t * math.cos(math.radians(x))        # (B6) Normal component of the absolute tip velocity [m/s]
Cm1rmsn = Cm1rms * math.cos(math.radians(x))    # (B6) Normal component of the root mean square velocity [m/s]
rho1h = rho0 * (1 - (Cm1h ** 2 / (2 * Cp * T00))) ** (1 / (k - 1))      # Inlet density of the air at the hub [kg/m^3]
rho1rms = rho0 * (1 - (Cm1rms ** 2 / (2 * Cp * T00))) ** (1 / (k - 1))  # Root mean square density of the air [kg/m^2]
rho1t = rho0 * (1 - (Cm1t ** 2 / (2 * Cp * T00))) ** (1 / (k - 1))      # Inlet density of the air at the blade tip [kg/m^3]
mdoto = 2 * math.pi * (((h0 + h1) / 6) * ((2 - (h1 / h0)) * (rho1h * rh1 * Cm1hn) + ((h0 + h1) ** 2 / (h0 * h1)) * (rho1rms * r1rms * Cm1rmsn) + (2 - (h0 / h1)) * (rho1t * rt1 * Cm1tn)))      # (B7) Off Design point mass flow rate [kg/s]
U1to = math.pi * No * Ndes * 2 * rt1 / 60       # Off design blade tip velocity [m/s]
U1ho = math.pi * No * Ndes * 2 * rh1 / 60       # Off design blade hub velocity [m/s]
U1rmso = ((U1to ** 2 + U1ho ** 2) / 2) ** 0.5   # Off design rms velocity [m/s]
W1ho = (Cm1h ** 2 + U1ho ** 2) ** 0.5           # Off design hub relative velocity [m/s]
W1to = (Cm1t ** 2 + U1to ** 2) ** 0.5           # Off design tip relative velocity [m/s]
W1rmso = ((W1ho ** 2 + W1to ** 2) / 2) ** 0.5   # rms relative velocity [m/s]
beta1t = []
beta1h = []
beta1rms = []
T1 = T00i - (Cm1rms ** 2) / (2 * Cp)            # Inlet static temperature [K]
for i in range (0, len(VOVCR)):
    beta1t.append(math.degrees(math.atan(Cm1t[i] / U1to)))              # Hub inlet relative angle [degrees]
    beta1h.append(math.degrees(math.atan(Cm1h[i] / U1ho)))              # Tip inlet relative angle [degrees]
    beta1rms.append(((beta1t[i] ** 2 + beta1h[i] ** 2) / 2) ** 0.5)     # rms inlet relative angle [degrees]


### Inducer Incidence Loss--------------------------------------------------------------------------------
BBF = 1 - (ZB * tu) / (2 * math.pi * r1rms)     # (41) Blade Blockage Factor [-]
eps = []
betaopt = []
WL = []
dhinc = []
T1orel = []
WCR = []
W1rmseff = []
T0T1 = []
T1a = []
T1rmso = []
P1rmso = []
P1arms = []
for i in range(0, len(VOVCR)):
    eps.append(math.degrees(math.atan((1 - BBF) * math.tan(math.radians(beta1rms[i])) / (1 + BBF * math.tan(math.radians(beta1rms[i])) ** 2))))     # (B40) Difference between compressor inlet relative flow angle and optimum incidence angle [degrees]
    betaopt.append(beta1rms[i] - eps[i])                    # (B42) Optimum relative flow angle [degrees]
    WL.append(W1rmso[i] * math.sin(math.radians(abs(betaopt[i] - beta1rms[i]))))    # (B43) Component of relative velocity lost [m/s]
    dhinc.append((WL[i] ** 2) / (2 * Cp))                   # (B44) Enthalpy loss due to incidence [J/kg]
    T1orel.append(T1[i] + W1rmso[i] ** 2 / (2 * Cp))        # Off design inlet relative temperature [K]
    WCR.append((2 * (k - 1)/(k + 1) * R * T1orel[i]) ** 0.5)            # Critical inlet relative velocity [m/s]
    W1rmseff.append(W1rmso[i] * math.cos(betaopt[i] - beta1rms[i]))     # Effective relative velocity [m/s]
    T0T1.append(1 - (k - 1) / (k + 1) * (W1rmseff[i] / WCR[i]) ** 2)    # Ratio of inlet static temperatures [-]
    T1a.append(T1orel[i] * T0T1[i])                         # Temperature just inside the blade [K]
    T1rmso.append(T00i[i] - (Cm1rms[i] ** 2 / (2 * Cp)))    # Off design Root mean square of static temperature at inlet [K]
    P1rmso.append(P00 * (T1rmso[i] / T00i[i]) ** (k / (k - 1)))             # Off design root mean square of static pressure at the inlet [Pa]
    P1arms.append(P1rmso[i] * math.exp(( - 1 * dhinc[i]) / (T1a[i] * R)))   # (B45) Total pressure just inside the bladed row [Pa]


### Impeller Work and Losses----------------------------------------------------------------------------
U2o = U1to / (rt1 / (D2 / 2))       # Off design exit blade velocity [m/s]
dhest = U2o ** 2                    # (B46) Initial approximation of enthalpy rise in impeller [J/kg]
T2oestabs = (dhest / (Cp * T00) + 1) * T00          # (B47) Estimate of the off design impeller exit total temperature [K]
rho2o = rho1 * (T2oestabs / T00) ** (1 / (k - 1))   # (B48) Off design impeller exit density [kg/m^3]


def Densityiteration(rho2o):
    """
    The Densityiteration(rho2o) function takes an initial guess of the impeller outlet density using
    Equation B48 above. It then uses this initial guess to calculate a series of velocities,
    temperatures and enthalpies corresponding to this initial guess before re-calculating the
    density.
    """
    Vm2m = mdoto / (math.pi * rho2o * D2 * b2)      # (B49) Meridional component of exit absolute velocity [m/s]
    VSL = U2o * (1 - sigma)                         # (B51) Slip velocity [m/s]
    Vtheta2 = (U2o - Vm2m * math.tan(math.radians( - beta2b)) - VSL)    # (B50) Tangential component of exit absolute velocity [m/s]
    T1orelrms = T1 + W1rmso ** 2 / (2 * Cp)         # Relative root mean square temperature [K]
    T2orel = T1orelrms + ((U2o ** 2 - U1t ** 2) / (2 * Cp))   # (B52) Exit temperature in the relative reference frame [K]
    #T2orel = T1 + ((U2o ** 2 - U1t ** 2) / (2 * Cp))
    Wtheta2 = U2o - Vtheta2                         # (B53) Tangential component of relative exit velocity [m/s]
    W2 = ((Vm2m ** 2) + (Wtheta2 ** 2)) ** 0.5      # (B54) Relative exit velocity [m/s]
    T2o = (T2orel - ((W2 ** 2) / (2 * Cp)))         # (B55) Off design point exit temperature [K]
    V2 = ((Vm2m ** 2) + (Vtheta2 ** 2)) ** 0.5      # (B56) Off design point absolute exit velocity [m/s]
    T2oabs = (T2o + (V2 ** 2) / (2 * Cp))           # (B57) Off design point exit temperature in the absolute reference frame [K]
    dhaero = (Cp * T00 * (T2oabs / T00 - 1))        # (B61) Aerodynamic enthalpy rise [J/kg]
    qaero = dhaero / (U2o ** 2)                     # (B60) Dimensionless actual head [-]
    Df = (1 - W2 / W1to + (kBL * qaero) / ((W1to / U2) * ((ZB / math.pi) * (1 - 2 * rt1 / D2) + 2 * 2 * rt1 / D2)))     # (B59)Diffusion factor [-]
    dhBL = (0.05 * Df ** 2 * U2o ** 2)              # (B58) Work loss due to blade loading [J/kg]
    Re = U2o * D2 * rho1rms / visc                  # (B63) Reynolds number of the exit flow [-]
    dhDF = (0.01356 * rho2o * U2o ** 3 * D2 ** 2 / (mdoto * Re ** 0.2))         # (B62) Impeller disk friction loss [J/kg]
    D1rms = math.sqrt(((( 2 * rt1) ** 2) + ((2 * rh1) ** 2)) / 2)               # Rootmean square of the diameter [m]
    Lendia = 0.5 * (1 - (D1rms / 0.3048)) / (math.cos(math.radians(beta2b)))    # (B65) Blade length to diameter ratio [-]
    HYDdia = 1 / (ZB / (math.pi * math.cos(math.radians(beta2b)) + D2 / b2)) + (2  * rt1 / D2) / (2 / (1 - k) + 2 * ZB / (math.pi * (1 + k)) * math.sqrt(1 + (math.tan(math.radians(beta1) **2 ) * (1 + k ** 2 / 2))))  # Ratio of hydraulic diameter and exit diameter [-]
    WRelExt = 0.5 * ((Cm1rms / U2o) ** 2 + (D1rms / D2) ** 2 + (W2 / W1to) ** 2 * ((Cm1rms / U2o) ** 2 + (2 * rt1 / D2) **2 ))      # (B67) Ratio of mean relative velocity and impeller exit velocity^2 [-]
    dhSF = ((kSF * Cf * Lendia * WRelExt * U2o ** 2) / HYDdia)      # (B64) Skin Friction loss [J/kg]
    dhid = (dhaero - dhinc - dhSF - dhDF - dhBL)    # (B68) Ideal enthalpy rise [J/kg]
    etaR = dhid / dhaero                            # (B69) Impeller efficiency [-]
    P2oabs = (P1arms * (etaR * dhaero / (Cp * T00) + 1) ** (k / (k - 1)))       # (B70) Iteration of the off design exit absolute pressure [Pa]
    for Temp in T2o:
        if Temp < 0:             # MSG: Added this if statement to raise error if negative temperatures
            raise ValueError("Temperature T2o is negative") 
    P2o = (P2oabs / ((T2oabs / T2o) ** (k / (k - 1))))      # (B71) Iteration of the off design exit pressure [Pa]      
    rho2oit = P2o / (R * T2o)                       # (B72) Iteration of the off design exit density [kg/m^3]
    return [rho2oit, T2o, dhaero, dhBL, dhDF, dhSF, dhid, T2oabs, P2oabs, Vtheta2, Vm2m, Df, P2o]


def Density():
    """
    The Density() function selects the value corresponding to the variables listed below when the
    output density from the DensityIteration(rho2o) function is within 0.1% of the input
    """
    rhoinit =rho2o
    RHO = [] 
    T2o = []
    dhaero = []
    dhBL = []
    dhDF = []
    dhSF = []
    dhid = []
    T2oabs = []
    P2oabs = []
    Vtheta2 = []
    Vm2m = []
    Df = []
    P2o = []
    for i in range(0, len(VOVCR)):
        rhoafter = Densityiteration(rho2o)[0][i]
        rho = [rhoinit, rhoafter]
        while abs((rho[- 1]) - (rho[ - 2])) > 0.00001:
            if abs((rho[- 1]) - (rho[- 2])) > 0.001:
                rho.append(Densityiteration(rho[- 1])[0][i])
            else:
                if len(RHO) < i + 1:
                    rho.append(Densityiteration(rho[- 1])[0][i])
                    RHO.append(Densityiteration(rho[- 1])[0][i])
                    T2o.append(Densityiteration(rho[- 1])[1][i])
                    dhaero.append(Densityiteration(rho[- 1])[2][i])
                    dhBL.append(Densityiteration(rho[- 1])[3][i])
                    dhDF.append(Densityiteration(rho[- 1])[4][i])
                    dhSF.append(Densityiteration(rho[- 1])[5][i])
                    dhid.append(Densityiteration(rho[- 1])[6][i] )
                    T2oabs.append(Densityiteration(rho[- 1])[7][i])
                    P2oabs.append(Densityiteration(rho[- 1])[8][i])
                    Vtheta2.append(Densityiteration(rho[- 1])[9][i])
                    Vm2m.append(Densityiteration(rho[- 1])[10][i])
                    Df.append(Densityiteration(rho[- 1])[11][i])
                    P2o.append(Densityiteration(rho[- 1])[12][i])
                else:
                    rho.append(Densityiteration(rho[- 1])[0][i])
    return [RHO, T2o, dhaero, dhBL, dhDF, dhSF, dhid, T2oabs, P2oabs, Vtheta2, Vm2m, Df, P2o]

rho2 = np.array(Density()[0])
T2o = np.array(Density()[1])
dhaero = np.array(Density()[2])
dhBL = np.array(Density()[3])
dhDF = np.array(Density()[4])
dhSF = np.array(Density()[5])
dhid = np.array(Density()[6])
T2oabs = np.array(Density()[7])
P2oabs = np.array(Density()[8])
Vtheta2 = np.array(Density()[9])
Vm2m = np.array(Density()[10])
Df = np.array(Density()[11])
P2o = np.array(Density()[12])
C2o = []

for i in range(0, len(VOVCR)):
    C2o.append((Vm2m[i] ** 2 + Vtheta2[i] ** 2) ** 0.5)


### Recirculation Loss-------------------------------------------------------------------------------------
dhRC = []
alpha2 = []
for i in range(0, len(VOVCR)):
    alpha2.append(math.degrees(math.atan(Vtheta2[i] / Vm2m[i])) )                               # Exit velocity flow angle [degrees]
    dhRC.append(0.02 * math.sqrt(math.tan(math.radians(alpha2[i]))) * Df[i] ** 2 * U2o ** 2)    # Enthalpy loss from recirculation [kJ/kg]


### Exit losses--------------------------------------------------------------------------------------------
etad = 0.85
CpDi = 1 - (AR ** - 2)
CpD = etad * CpDi
M3o = []
P3oabs = []
dhVLD = []
P3o = []
C3o = []
P03o = []
for i in range(0, len(VOVCR)):
    P3o.append(CpD * 0.5 * rho2[i] * (C2o[i]) ** 2 + P2o[i])
    C3o.append(C2o[i] / AR)
    P03o.append(P3o[i] + 0.5 * rho2[i] * C3o[i] ** 2)
    M3o.append(C3o[i] / (math.sqrt(k * R * T2oabs[i])))
    P3oabs.append(P3o[i] * (1 + (k - 1) / 2 * M3o[i] ** 2) ** (k / (k - 1)))    # (B85) Absolute diffuser throat pressure [Pa]
    dhVLD.append(Cp * T2oabs[i] * ((P3o[i] / P3oabs[i]) ** ((k - 1) / k) - (P3o[i] / P2oabs[i]) ** ((k - 1) / k)))  # (B86) Vaneless diffuser loss [kJ/kg]


### Overall performance-----------------------------------------------------------------------------------
etao = []
Pro = []
for i in range(0, len(VOVCR)):
    etao.append((dhaero[i] - (dhinc[i] + dhBL[i] + dhSF[i] + dhVLD[i])) / (dhaero[i] + dhRC[i] + dhDF[i]))
    Pro.append(P3o[i] / P00)

print("Mass flow rate =", mdoto)
print("Pressure ratio =", Pro)
print("Efficiency =", etao)

plt.figure()
plt.plot(mdoto, etao, 'b-')
plt.xlabel("Mass flow rate (kg/s)")
plt.ylabel("Efficiency (%)")
plt.title("Efficiency contour")
plt.xlim(0, 1)
plt.show()

plt.figure()
plt.plot(mdoto, Pro, 'r', label = "120000 RPM", linewidth = 1.25)
plt.xlabel("Mass flow rate (kg/s)")
plt.ylabel("Pressure ratio (-)")
plt.title("Compressor Map")


### Matching-------------------------------------------------------------------------------------------
comp_power = [100]      # Compressor power [kW]
speed = [120]           # Shaft Speed [krpm]
T1 = 293                # Inlet stagnation temperature [K]
cp = 1.005              # Specific heat at constant pressure [kJ]
y = 1.4                 # Ratio of specific heats
P = 0.1013              # Inlet pressure [MPa]
eta = 0.6               # Assumed efficiency


### Line of constant power------------------------------------------------------------------
def pressure_ratio1(mdot, eta, comp_power, T1, y, cp):
    return (eta * comp_power / (mdot * cp * T1) + 1) ** (y / (y - 1))

mdot = np.linspace(0.0001, 1.4, 60)
for i in range(1):
    Pr = pressure_ratio1(mdot, eta, comp_power[i], T1, y, cp)
plt.plot(mdot, Pr, 'g-', label = "100 kW", linewidth = 1.25)
plt.xlim(0, 0.6)
plt.ylim(1, 8)
plt.legend()
plt.show()