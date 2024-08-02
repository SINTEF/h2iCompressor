"""
The following Python code uses the inducer and impeller geometryV8 from the geometryV8 code
(geometryV8.py) and predicts the performance of the compressor at off-design points.

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
import settingsOffDesign
# import logic
import geometryV8     # MSG: This runs the geometryV8.py file. Avoid this by including main function?
import math
import numpy as np
import matplotlib.pyplot as plt

# Set plot parameters------------------------------------------------------------------------------
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams.update({'font.size': 15})
plt.close('all')

print('Running off-design_performance.py')

def off_design_performance(Nrpm):
    # Control, Determining point of interest -----------------
    bladeAngle = settingsOffDesign.bladeAngle
    bladeNumber = settingsOffDesign.bladeNumber

    iB = int(np.where(geometryV8.beta2BArr == bladeAngle)[0])
    iZ = int(np.where(geometryV8.ZBarr == bladeNumber)[0])
    beta2B = geometryV8.beta2BArr[iB]
    b2 = geometryV8.b2Mat[iB][iZ]
    ZB = geometryV8.ZBarr[iZ]
    sigma = geometryV8.sigmaMat[iB][iZ]

    ### Parameters, non-variable ------------------------------------------------------------------------------------------------
    # Ndes = geometryV8.Ndes
    # Ndes = logic.N
    Ndes = Nrpm
    # Ndes = logic.N
    No = 1.0                                                # Off design speed percentage that requires performance prediction [rpm]
    tu = 0.002                                              # Blade thickness [m]
    # rho0 = 1.1839                                           # Stagnation density of air [kg/m^3]
    rho0 = settingsOffDesign.rho0                           # Stagnation density of H2 [kg/m^3]
    curvet1 = 0                                             # Inducer inlet tip wall curvature [m^-1]
    curveh1 = 0                                             # Inducer inlet hub wall curvature [m^-1]
    x = 10                                                  # Streamline angle [degrees] from axial direction
    VelCr = math.sqrt(2 * settingsOffDesign.k / (settingsOffDesign.k + 1) * settingsOffDesign.R * settingsOffDesign.T00)      # Critical Velocity wrt resonance [m/s]  
    V0DivVcr = np.linspace(0.1, 0.686, 50)                     # Compressor inlet absolute critical velocity ratio wrt resonance [-]
    Cm1h = V0DivVcr * VelCr                                    # Absolute meridional velocity of the hub [m/s]
    kBL = 0.6                                               # Blading loss coefficient [-]
    visc = 1.778e-5                                         # Air viscosity based on total conditions
    visc = 0.88e-5                                          # H2 viscosity based on total conditions
    kSF = 7.0                                               # Skin friction coefficient [-]
    Cf = 0.01                                               # Friction coefficient [-]
    T00i = np.full(len(V0DivVcr), settingsOffDesign.T00)


    # ---------------------------- Calculation of Inlet Velocity Triangles and Compressor Weight Flow (Swirl free) ----------------------------

    curve1rms = math.sqrt((curvet1 ** 2 + curveh1 ** 2) / 2)        # Root mean square of the inducer inlet hub and tip wall curvature [m^-1]
    r1rms = math.sqrt((geometryV8.r1 ** 2 + geometryV8.rh1 ** 2) / 2)  # Root mean square of the hub and tip inlet radius [m]

    # Geometrical inlet properties
    h0 = r1rms - geometryV8.rh1            # (B4) Spacing for numerical integration [-]
    h1 = geometryV8.r1 - r1rms            # (B5) Spacing for numerical integration [-]

    # Absolute meridonal velocities
    Cm1rms = Cm1h * math.exp((h0 / 2) * (curveh1 + curve1rms))      # (B2) Absolute meridional root mean square velocity [m/s]
    Cm1t = Cm1h * math.exp(((h0 + h1) / 6) * ((2 - (h1 / h0)) * curveh1 + (((h0 + h1) ** 2) / (h1 * h0)) * curve1rms + (2 - (h0 / h1)) * curvet1))      # (B3) Absolute meridional velocity of the tip [m/s]

    # Velocities through the normal flow area at inlet
    Cm1hn = Cm1h * math.cos(math.radians(x))        # (B6) Normal component of the absolute hub velocity [m/s], through flow area
    Cm1tn = Cm1t * math.cos(math.radians(x))        # (B6) Normal component of the absolute tip velocity [m/s], through flow area
    Cm1rmsn = Cm1rms * math.cos(math.radians(x))    # (B6) Normal component of the root mean square velocity [m/s], through flow area

    # Densities from stagnation temperature and different velocites
    rho1h = rho0 * (1 - (Cm1h ** 2 / (2 * settingsOffDesign.Cp * settingsOffDesign.T00))) ** (1 / (settingsOffDesign.k - 1))      # Inlet density of the air at the hub [kg/m^3], from stagnation temperature
    rho1rms = rho0 * (1 - (Cm1rms ** 2 / (2 * settingsOffDesign.Cp * settingsOffDesign.T00))) ** (1 / (settingsOffDesign.k - 1))  # Root mean square density of the air [kg/m^2], 
    rho1t = rho0 * (1 - (Cm1t ** 2 / (2 * settingsOffDesign.Cp * settingsOffDesign.T00))) ** (1 / (settingsOffDesign.k - 1))      # Inlet density of the air at the blade tip [kg/m^3]

    # mass flow rate
    mdoto = 2 * math.pi * (((h0 + h1) / 6) * ((2 - (h1 / h0)) * (rho1h * geometryV8.rh1 * Cm1hn) + ((h0 + h1) ** 2 / (h0 * h1)) * (rho1rms * r1rms * Cm1rmsn) + (2 - (h0 / h1)) * (rho1t * geometryV8.r1 * Cm1tn)))      # (B7) Off Design point mass flow rate [kg/s]

    # Blade speeds
    U1to = math.pi * No * Ndes * 2 * geometryV8.r1 / 60       # Off design blade tip velocity [m/s]
    U1ho = math.pi * No * Ndes * 2 * geometryV8.rh1 / 60       # Off design blade hub velocity [m/s]
    U1rmso = ((U1to ** 2 + U1ho ** 2) / 2) ** 0.5   # Off design rms velocity [m/s]

    # Relative velocities at inlet
    W1ho = (Cm1h ** 2 + U1ho ** 2) ** 0.5           # Off design hub relative velocity [m/s]
    W1to = (Cm1t ** 2 + U1to ** 2) ** 0.5           # Off design tip relative velocity [m/s]
    W1rmso = ((W1ho ** 2 + W1to ** 2) / 2) ** 0.5   # rms relative velocity [m/s]

    # Relative flow inlet angles beta (not blade angle, flow angle!)
    beta1t = []             # Hub inlet relative angle [degrees]
    beta1h = []             # Tip inlet relative angle [degrees]
    beta1rms = []           # rms inlet relative angle [degrees]

    T1 = T00i - (Cm1rms ** 2) / (2 * settingsOffDesign.Cp)   # Inlet static temperature [K]
    for i in range (0, len(V0DivVcr)):
        beta1t.append(math.degrees(math.atan(Cm1t[i] / U1to)))              # Hub inlet relative angle [degrees]
        beta1h.append(math.degrees(math.atan(Cm1h[i] / U1ho)))              # Tip inlet relative angle [degrees]
        beta1rms.append(((beta1t[i] ** 2 + beta1h[i] ** 2) / 2) ** 0.5)     # rms inlet relative angle [degrees]



    ###---------------------------------------- Inducer Incidence Loss----------------------------------------
    # declaring variables before use
    BBF = 1 - (ZB * tu) / (2 * math.pi * r1rms)     # (41) Blade Blockage Factor [-]        
    eps = []                        # Difference between compressor inlet relative flow angle and optimum incidence angle [degrees]
    betaOpt = []                    # Optimum relative flow angle [degrees]
    WL = []                         # Relative velocity loss [m/s]
    dhInc = []                      # Enthalpy loss due to incidence [J/kg]
    T1oRel = []                     # Off design inlet relative temperature [K]
    W1cr = []                        # Critical inlet relative velocity [m/s]
    W1rmsEff = []                   # Effective relative velocity [m/s]
    T0divT1 = []                       # Ratio of inlet static temperatures [-]
    T1a = []                        # Temperature just inside the blade [K]
    T1rmso = []                     # Off design Root mean square of static temperature at inlet [K]
    P1rmso = []                     # Off design root mean square of static pressure at the inlet [Pa]
    P1arms = []                     # Total pressure just inside the bladed row [Pa]

    # Looping to find every variable of interest
    for i in range(0, len(V0DivVcr)):
        eps.append(math.degrees(math.atan((1 - BBF) * math.tan(math.radians(beta1rms[i])) / (1 + BBF * math.tan(math.radians(beta1rms[i])) ** 2))))     # (B40) Difference between compressor inlet relative flow angle and optimum incidence angle [degrees]
        betaOpt.append(beta1rms[i] - eps[i])                                                                                    # (B42) Optimum relative flow angle [degrees]
        WL.append(W1rmso[i] * math.sin(math.radians(abs(betaOpt[i] - beta1rms[i]))))                                            # (B43) Component of relative velocity lost [m/s]
        
        # THIS ONE IS FISHY WRT UNITS OF J/KG, UNECCESARY CP?:
        # dhInc.append((WL[i] ** 2) / (2 * settingsOffDesign.Cp))                                                                 # (B44) Enthalpy loss due to incidence [J/kg] => THIS EQUATION GIVES KELVIN AS UNITS, WHICH IS NOT ALIGNED WITH FURTHER CALCULATIONS
        dhInc.append((WL[i] ** 2) / (2 )) #* settingsOffDesign.Cp))                                                                 # (B44) Enthalpy loss due to incidence [J/kg] => THIS EQUATION GIVES KELVIN AS UNITS, WHICH IS NOT ALIGNED WITH FURTHER CALCULATIONS
        
        T1oRel.append(T1[i] + W1rmso[i] ** 2 / (2 * settingsOffDesign.Cp))                                                      # Off design inlet relative temperature [K]
        W1cr.append((2 * (settingsOffDesign.k - 1)/(settingsOffDesign.k + 1) * settingsOffDesign.R * T1oRel[i]) ** 0.5)          # Critical inlet relative velocity [m/s]
        W1rmsEff.append(W1rmso[i] * math.cos(betaOpt[i] - beta1rms[i]))                                                         # Effective relative velocity [m/s]
        T0divT1.append(1 - (settingsOffDesign.k - 1) / (settingsOffDesign.k + 1) * (W1rmsEff[i] / W1cr[i]) ** 2)                 # Ratio of inlet static temperatures [-]
        T1a.append(T1oRel[i] * T0divT1[i])                                                                                      # Temperature just inside the blade [K]
        T1rmso.append(T00i[i] - (Cm1rms[i] ** 2 / (2 * settingsOffDesign.Cp)))                                                  # Off design Root mean square of static temperature at inlet [K], STAGNATION EQUATION 
        P1rmso.append(settingsOffDesign.P00 * (T1rmso[i] / T00i[i]) ** (settingsOffDesign.k / (settingsOffDesign.k - 1)))       # Off design root mean square of static pressure at the inlet [Pa], ISENTROPIC PROPERTY
        P1arms.append(P1rmso[i] * math.exp(( - 1 * dhInc[i]) / (T1a[i] * settingsOffDesign.R)))                                 # (B45) Total pressure just inside the bladed row [Pa]

    # print('debug')

    ### Impeller Work and Losses----------------------------------------------------------------------------
    U2o = U1to / (geometryV8.r1 / (geometryV8.r2))                                                                      # Off design exit blade velocity [m/s], CONSTANT ANGULAR VELOCITY
    dhest = U2o ** 2                                                                                                    # (B46) Initial approximation of enthalpy rise in impeller [J/kg], EULER EQUATION APPROXIMATION zero inlet swirl
    T2oEstAbs = (dhest / (settingsOffDesign.Cp * settingsOffDesign.T00) + 1) * settingsOffDesign.T00                    # (B47) Estimate of the off design impeller exit total temperature [K], ENERGY/EULER EQUATION
    rho2o = geometryV8.rho1 * (T2oEstAbs / settingsOffDesign.T00) ** (1 / (settingsOffDesign.k - 1))                    # (B48) Off design impeller exit density [kg/m^3], ISENTROPIC RELATION


    def Densityiteration(rho2o):
        """
        The Densityiteration(rho2o) function takes an initial guess of the impeller outlet density using
        Equation B48 above. It then uses this initial guess to calculate a series of velocities,
        temperatures and enthalpies corresponding to this initial guess before re-calculating the
        density.
        """
        Vm2m = mdoto / (math.pi * rho2o * (2*geometryV8.r2) * b2)                        # (B49) Meridional component of exit absolute velocity [m/s]  , MASS BALANCE                
        VSL = U2o * (1 - sigma)                                                                     # (B51) Slip velocity [m/s]  , SLIP RELATION          
        Vtheta2 = (U2o - Vm2m * math.tan(math.radians( - beta2B)) - VSL)                            # (B50) Tangential component of exit absolute velocity [m/s]    , VELOCITY TRIANGLE       
        T1orelrms = T1 + W1rmso ** 2 / (2 * settingsOffDesign.Cp)                                   # Relative root mean square temperature [K]   , STAGNATION TEMPERATURE
        T2orel = T1orelrms + ((U2o ** 2 - geometryV8.U1t ** 2) / (2 * settingsOffDesign.Cp))        # (B52) Exit temperature in the relative reference frame [K]  , STAGNATION RELATION
        T2orel = T1orelrms + ((U2o ** 2 - U1to ** 2) / (2 * settingsOffDesign.Cp))        # (B52) Exit temperature in the relative reference frame [K]  , STAGNATION RELATION
        #T2orel = T1 + ((U2o ** 2 - U1t ** 2) / (2 * Cp))
        Wtheta2 = U2o - Vtheta2                                                                     # (B53) Tangential component of relative exit velocity [m/s], VELOCITY TRIANGLE
        W2 = ((Vm2m ** 2) + (Wtheta2 ** 2)) ** 0.5                                                  # (B54) Relative exit velocity [m/s]   ,  VELOCITY TRIANGLE
        T2o = (T2orel - ((W2 ** 2) / (2 * settingsOffDesign.Cp)))                                   # (B55) Off design point exit temperature [K]  , STAGNATION RELATION
        V2 = ((Vm2m ** 2) + (Vtheta2 ** 2)) ** 0.5                                                  # (B56) Off design point absolute exit velocity [m/s]  ,  VELOCITY TRIANGLE  
        T2oabs = (T2o + (V2 ** 2) / (2 * settingsOffDesign.Cp))                                     # (B57) Off design point exit temperature in the absolute reference frame [K]   , STAGNATION RELATION
        dhaero = (settingsOffDesign.Cp * settingsOffDesign.T00 * (T2oabs / settingsOffDesign.T00 - 1))        # (B61) Aerodynamic enthalpy rise [J/kg]  , Cp*dT
        qaero = dhaero / (U2o ** 2)                                                                           # (B60) Dimensionless actual head [-]   , WORK COEFFICIENT
        Df = (1 - W2 / W1to + (kBL * qaero) / ((W1to / geometryV8.U2) * ((ZB / math.pi) * (1 - 2 * geometryV8.r1 / geometryV8.D2) + 2 * 2 * geometryV8.r1 / geometryV8.D2)))     # (B59)Diffusion factor [-]  , EMPIRICAL
        dhBL = (0.05 * Df ** 2 * U2o ** 2)                                                              # (B58) Work loss due to blade loading [J/kg]   , ----
        Re = U2o * geometryV8.D2 * rho1rms / visc                                                       # (B63) Reynolds number of the exit flow [-]   ,  ----
        dhDF = (0.01356 * rho2o * U2o ** 3 * geometryV8.D2 ** 2 / (mdoto * Re ** 0.2))                  # (B62) Impeller disk friction loss [J/kg]   , ----
        D1rms = math.sqrt(((( 2 * geometryV8.r1) ** 2) + ((2 * geometryV8.rh1) ** 2)) / 2)              # Rootmean square of the diameter [m]  , RMS
        LenDivDia = 0.5 * (1 - (D1rms / 0.3048)) / (math.cos(math.radians(beta2B)))                    # (B65) Blade length to diameter ratio [-]   ,  ----
        HydDiaDivExitDia = 1 / (ZB / (math.pi * math.cos(math.radians(beta2B)) + geometryV8.D2 / b2)) + (2  * geometryV8.r1 / geometryV8.D2) / (2 / (1 - settingsOffDesign.k) + 2 * ZB / (math.pi * (1 + settingsOffDesign.k)) * math.sqrt(1 + (math.tan(math.radians(geometryV8.beta1) **2 ) * (1 + settingsOffDesign.k ** 2 / 2))))  # Ratio of hydraulic diameter and exit diameter [-]
        WRelDivWext = 0.5 * ((Cm1rms / U2o) ** 2 + (D1rms / geometryV8.D2) ** 2 + (W2 / W1to) ** 2 * ((Cm1rms / U2o) ** 2 + (2 * geometryV8.r1 / geometryV8.D2) **2 ))      # (B67) Ratio of mean relative velocity and impeller exit velocity^2 [-]
        dhSF = ((kSF * Cf * LenDivDia * WRelDivWext * U2o ** 2) / HydDiaDivExitDia)             # (B64) Skin Friction loss [J/kg]
        dhid = (dhaero - dhInc - dhSF - dhDF - dhBL)                                            # (B68) Ideal enthalpy rise [J/kg]
        etaR = dhid / dhaero                                                                    # (B69) Impeller efficiency [-]
        P2oabs = (P1arms * (etaR * dhaero / (settingsOffDesign.Cp * settingsOffDesign.T00) + 1) ** (settingsOffDesign.k / (settingsOffDesign.k - 1)))       # (B70) Iteration of the off design exit absolute pressure [Pa]
        for Temp in T2o:
            if Temp < 0:             # MSG: Added this if statement to raise error if negative temperatures
                raise ValueError("Temperature T2o is negative", Temp) 
            
        P2o = (P2oabs / ((T2oabs / T2o) ** (settingsOffDesign.k / (settingsOffDesign.k - 1))))          # (B71) Iteration of the off design exit pressure [Pa]      
        rho2oit = P2o / (settingsOffDesign.R * T2o)                                                     # (B72) Iteration of the off design exit density [kg/m^3]
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
        for i in range(0, len(V0DivVcr)):
            rhoafter = Densityiteration(rho2o)[0][i]
            rho = [rhoinit, rhoafter]                               # Control array for iteration
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

    # Enthalpy loss and power contributions 
    dhaero = np.array(Density()[2])             # Aerodynamic enthalpy rise
    dhBL = np.array(Density()[3])               # Blade loading
    dhDF = np.array(Density()[4])               # Disk friction
    dhSF = np.array(Density()[5])               # Skin friciton
    dhid = np.array(Density()[6])               # Ideal enthaply riise

    T2oabs = np.array(Density()[7])         # Stagnation temperature at exit
    P2oabs = np.array(Density()[8])         # Stagnation pressure at exit
    Vtheta2 = np.array(Density()[9])        # Tangential velocity component at exit
    Vm2m = np.array(Density()[10])          # Meridonal velocity component at exit
    Df = np.array(Density()[11])            # Diffusion factor
    P2o = np.array(Density()[12])           # Stagnation pressure at exit after iteration
    C2o = []                                # Absolute flow velocity at oulet

    for i in range(0, len(V0DivVcr)):
        C2o.append((Vm2m[i] ** 2 + Vtheta2[i] ** 2) ** 0.5)


    ### Recirculation Loss, work done on fluid going back --------------------------------------------------------------------------
    dhRC = []
    alpha2 = []
    for i in range(0, len(V0DivVcr)):
        alpha2.append(math.degrees(math.atan(Vtheta2[i] / Vm2m[i])) )                               # Exit velocity flow angle [degrees]
        dhRC.append(0.02 * math.sqrt(math.tan(math.radians(alpha2[i]))) * Df[i] ** 2 * U2o ** 2)    # Enthalpy loss from recirculation [kJ/kg]


    ### Exit losses--------------------------------------------------------------------------------------------
    etad = settingsOffDesign.etad
    CpDi = 1 - (settingsOffDesign.AR ** - 2)
    CpD = etad * CpDi
    M3o = []
    P3oabs = []
    dhVLD = []
    P3o = []
    C3o = []
    P03o = []
    for i in range(0, len(V0DivVcr)):
        P3o.append(CpD * 0.5 * rho2[i] * (C2o[i]) ** 2 + P2o[i])
        C3o.append(C2o[i] / settingsOffDesign.AR)
        P03o.append(P3o[i] + 0.5 * rho2[i] * C3o[i] ** 2)
        M3o.append(C3o[i] / (math.sqrt(settingsOffDesign.k * settingsOffDesign.R * T2oabs[i])))
        P3oabs.append(P3o[i] * (1 + (settingsOffDesign.k - 1) / 2 * M3o[i] ** 2) ** (settingsOffDesign.k / (settingsOffDesign.k - 1)))    # (B85) Absolute diffuser throat pressure [Pa]
        dhVLD.append(settingsOffDesign.Cp * T2oabs[i] * ((P3o[i] / P3oabs[i]) ** ((settingsOffDesign.k - 1) / settingsOffDesign.k) - (P3o[i] / P2oabs[i]) ** ((settingsOffDesign.k - 1) / settingsOffDesign.k)))  # (B86) Vaneless diffuser loss [kJ/kg]


    ### Overall performance-----------------------------------------------------------------------------------
    etao = []
    Pro = []
    for i in range(0, len(V0DivVcr)):
        etaAppend = (dhaero[i] - (dhInc[i] + dhBL[i] + dhSF[i] + dhVLD[i])) / (dhaero[i] + dhRC[i] + dhDF[i])       #debugging
        etao.append(etaAppend)      
        P3append = P3o[i]           #debugging
        Pro.append(P3append / settingsOffDesign.P00)
        if P3append < settingsOffDesign.P00:
            debug= 0

    # print("\nMass flow rate =", mdoto)
    # print("\nPressure ratio =", Pro)
    # print("\nEfficiency =", etao)

    # plt.figure()
    # plt.plot(mdoto, etao, 'b-')
    # plt.xlabel("Mass flow rate (kg/s)")
    # plt.ylabel("Efficiency (%)")
    # plt.title("Efficiency contour")
    # # plt.xlim(0, 1)
    # plt.grid()

    # plt.figure()
    # plt.plot(mdoto, Pro, 'r', label = str(Ndes) +" RPM", linewidth = 2)
    # plt.xlabel("Mass flow rate (kg/s)")
    # plt.ylabel("Pressure ratio (-)")
    # plt.title("Compressor Map")


    ### Matching-------------------------------------------------------------------------------------------
    # comp_power = [100]      # Compressor power [kW]
    # comp_power = [round(geometryV8.WxMat[iB][iZ], 2) ]      # Compressor power [kW]
    # # speed = [120]           # Shaft Speed [krpm]    MSG: Not in use
    # T1 = settingsOffDesign.T00                # Inlet stagnation temperature [K]
    # cp = settingsOffDesign.Cp              # Specific heat at constant pressure [kJ]
    # y = settingsOffDesign.k                 # Ratio of specific heats
    # # P = 0.1013              # Inlet pressure [MPa]
    # eta = geometryV8.etaMat[iB][iZ]              # Assumed efficiency


    ### Line of constant power------------------------------------------------------------------------------
    # def pressure_ratio1(mdot, eta, comp_power, T1, y, cp):
    #     return (eta * comp_power / (mdot * cp * T1) + 1) ** (y / (y - 1))

    # mdot = np.linspace(0.0001, 1.4, 60)
    # for i in range(len(comp_power)):
    #     Pr = pressure_ratio1(mdot, eta, comp_power[i], T1, y, cp)
    #     plt.plot(mdot, Pr, 'g-', label = str(comp_power[i]) + " kW", linewidth = 2)
    # # plt.xlim(0, 0.6)
    # plt.ylim(1, 2)
    # plt.legend()
    # plt.grid()


    # print(str(geometryV8.r1) + '\n')
    # print(str(geometryV8.r2))

    print('Current RPM: ' + str(round(Ndes, 2)) + '     Inlet blade tip speed: ' + str(round(U1ho, 2)) + ' m/s     Outlet blade tip speed: ' + str(round(U2o, 2)) + ' m/s' )
    # print('Current RPM: ' + str(Ndes))

    # plt.show()


    return Pro, P03o, T2oabs, mdoto, etao