"""
The following Python code uses the inducer and impeller geometry obtainted from the geometry.py
and predicts the off-design stage performance for a given blade number and blade angle.

The method follows the one in Galvas (1973) with a few simplifications as described in Lowe (2016).

Note that the letter B and number inside the brackets in some of the lines of code correspond
to the equations in NASA’s Fortran model (B46, B47…B60 for example). The B denotes
Appendix B in the original paper and the number corresponds to the equation number in
approximate order of use in the code (Galvas, 1973). There is a jump in numbers between B7-
B40 and B72-B84. The missing equations correspond to inlet guide vane losses and vaned
diffuser losses, both of which were omitted for reasons outlined in Section 4.2 of Lowe (2016). 
Many of the variables and equations used in this code have not been defined in this thesis. 
The full set of equations and a description of the calculation process can be found in Appendix B of NASA’s
Fortran model.

References: Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244
            Galvas, Michael R. (1973). Fortran program for predicting off-design performance of centrifugal compressors, https://ntrs.nasa.gov/citations/19740001912 

TO-DO: 
        - Comments marked with "MSG" must be checked out 
        - Check that simplifications of Galvas (1973) made in Lowe (2016) make sense

Author(s): Petter Resell (summer intern, 2024), Martin Spillum Grønli (SINTEF Energy Research, 2025)
"""


# Import
import numpy as np
import matplotlib.pyplot as plt


def off_design_performance_rpm(Nrpm, Compressor, Fluid, InletConditions, IterationMatrix):
    """ Calculate stage off-design performance for a given rotational speed. To obtain the impeller outlet density, two functions are defined inside this function """
    
    print('Blade number =',Compressor.bladeNumber)                                  # Blade number for which off-design performance is calculated 
    print('Blade angle =', np.rad2deg(Compressor.bladeAngle))                       # Blade angle for which off-design performance is calculated

    iB = int(np.where(IterationMatrix.beta2B == Compressor.bladeAngle)[0])          # Index for blade angle
    iZ = int(np.where(IterationMatrix.ZB == Compressor.bladeNumber)[0])             # Index for blade number
    beta2B = IterationMatrix.beta2B[iB]                                             # Blade angle 
    ZB = IterationMatrix.ZB[iZ]                                                     # Blade number
    b2 = IterationMatrix.b2[iB][iZ]                                                 # Compressor outlet height 
    sigma = IterationMatrix.sigma[iB][iZ]                                           # Slip factor        

    if np.isnan(b2):
        raise Exception(f"Blade number and blade angle is not found to be among valid cases! \n Chose another point or adjust parameters!")

    # Inlet conditions
    VelCr = np.sqrt(2 * Fluid.k / (Fluid.k + 1) * Fluid.R * InletConditions.T00)    # Critical velocity [m/s]
    Cm1h = Compressor.V0DivVcr * VelCr                                              # Absolute meridional velocity of the hub [m/s]
    T00i = np.full(len(Compressor.V0DivVcr), InletConditions.T00)                   # Inlet stagnation temperature array [K]

    # Calculation of inlet velocity triangles and compressor weight flow (swirl free)
    curve1rms = np.sqrt((Compressor.curvet1 ** 2 + Compressor.curveh1 ** 2) / 2)    # Root mean square of the inducer inlet hub and tip wall curvature [m^-1]
    r1rms = np.sqrt((Compressor.r1 ** 2 + Compressor.rh ** 2) / 2)                  # Root mean square of the hub and tip inlet radius [m] 

    # Spacing for numerical integration
    h0 = r1rms - Compressor.rh                                                      # (B4) Spacing for numerical integration [-]
    h1 = Compressor.r1 - r1rms                                                      # (B5) Spacing for numerical integration [-]

    # Absolute meridonal velocities
    Cm1rms = Cm1h * np.exp((h0 / 2) * (Compressor.curveh1 + curve1rms))             # (B2) Absolute meridional root mean square velocity [m/s]
    Cm1t = Cm1h * np.exp(((h0 + h1) / 6) * ((2 - (h1 / h0)) * Compressor.curveh1 + (((h0 + h1) ** 2) / (h1 * h0)) * curve1rms + (2 - (h0 / h1)) * Compressor.curvet1))      # (B3) Absolute meridional velocity of the tip [m/s]

    # Velocities through the normal flow area at inlet
    Cm1hn = Cm1h * np.cos(np.radians(Compressor.x))                                 # (B6) Normal component of the absolute hub velocity through flow area [m/s] 
    Cm1tn = Cm1t * np.cos(np.radians(Compressor.x))                                 # (B6) Normal component of the absolute tip velocity through flow area [m/s]
    Cm1rmsn = Cm1rms * np.cos(np.radians(Compressor.x))                             # (B6) Normal component of the root mean square velocity through flow area [m/s]

    # Densities from stagnation temperature and different velocites
    rho1h = InletConditions.rho0 * (1 - (Cm1h ** 2 / (2 * Fluid.Cp * InletConditions.T00))) ** (1 / (Fluid.k - 1))                # Inlet density of the fluid at the hub [kg/m^3]
    rho1rms = InletConditions.rho0 * (1 - (Cm1rms ** 2 / (2 * Fluid.Cp * InletConditions.T00))) ** (1 / (Fluid.k - 1))            # Root mean square density of the fluid [kg/m^2] 
    rho1t = InletConditions.rho0 * (1 - (Cm1t ** 2 / (2 * Fluid.Cp * InletConditions.T00))) ** (1 / (Fluid.k - 1))                # Inlet density of the fluid at the blade tip [kg/m^3]

    mdoto = 2 * np.pi * (((h0 + h1) / 6) * ((2 - (h1 / h0)) * (rho1h * Compressor.rh * Cm1hn) + ((h0 + h1) ** 2 / (h0 * h1)) * (rho1rms * r1rms * Cm1rmsn) + (2 - (h0 / h1)) * (rho1t * Compressor.r1 * Cm1tn)))      # (B7) Off Design point mass flow rate [kg/s]
    
    # Blade speeds
    U1to = np.pi * Compressor.No * Nrpm * 2 * Compressor.r1 / 60                    # Off design blade tip velocity [m/s]
    U1ho = np.pi * Compressor.No * Nrpm * 2 * Compressor.rh / 60                    # Off design blade hub velocity [m/s]
    U1rmso = ((U1to ** 2 + U1ho ** 2) / 2) ** 0.5                                   # Off design root mean square velocity [m/s]

    # Relative velocities at inlet
    W1ho = (Cm1h ** 2 + U1ho ** 2) ** 0.5                                           # Off design hub relative velocity [m/s]  
    W1to = (Cm1t ** 2 + U1to ** 2) ** 0.5                                           # Off design tip relative velocity [m/s]  
    W1rmso = ((W1ho ** 2 + W1to ** 2) / 2) ** 0.5                                   # Root mean square relative velocity [m/s]

    # Relative flow inlet angles (not blade angle, but flow angle)
    beta1t = []                                                                     # Hub inlet relative angle [deg]
    beta1h = []                                                                     # Tip inlet relative angle [deg]
    beta1rms = []                                                                   # Root mean square inlet relative angle [deg]

    T1 = T00i - (Cm1rms ** 2) / (2 * Fluid.Cp)                                      # Inlet static temperature [K]
    for i in range (0, len(Compressor.V0DivVcr)):
        #beta1t.append(np.degrees(np.arctan(Cm1t[i] / U1to)))                       # Tip inlet relative angle [deg]    MSG: Should this be U1to/Cm1t instead? Seems to be different angle than what is used when finding eps. See https://ntrs.nasa.gov/api/citations/19930083715/downloads/19930083715.pdf. Also, is the tangential velocity of the fluid at the impeller tip inlet upstream of the blade assumed to be zero? Ref. eq. B26 and B25
        beta1t.append(np.degrees(np.arctan(U1to / Cm1t[i])))
        #beta1h.append(np.degrees(np.arctan(Cm1h[i] / U1ho)))                       # Hub inlet relative angle [deg]    MSG: Should this be U1to/Cm1t instead? Seems to be different angle than what is used when finding eps. See https://ntrs.nasa.gov/api/citations/19930083715/downloads/19930083715.pdf. Also, is the tangential velocity of the fluid at the impeller tip inlet upstream of the blade assumed to be zero? Ref. eq. B26 and B25
        beta1h.append(np.degrees(np.arctan(U1ho / Cm1h[i])))

        beta1rms.append(((beta1t[i] ** 2 + beta1h[i] ** 2) / 2) ** 0.5)             # Root mean square inlet relative angle [deg]

    # Inducer incidence loss
    BBF = 1 - (ZB * Compressor.tu) / (2 * np.pi * r1rms)                            # (41) Blade Blockage Factor [-]  
    print('Blade blockage factor (should be around 0.9 to align with Galvas, 1973):', BBF)       
    # Declaring variables before use
    eps = []                                                                        # Difference between compressor inlet relative flow angle and optimum incidence angle [deg]
    betaOpt = []                                                                    # Optimum relative flow angle [deg]
    WL = []                                                                         # Relative velocity loss [m/s]
    dhInc = []                                                                      # Enthalpy loss due to incidence [J/kg]
    T1oRel = []                                                                     # Off design inlet relative temperature [K]
    W1cr = []                                                                       # Critical inlet relative velocity [m/s]
    W1rmsEff = []                                                                   # Effective relative velocity [m/s]
    T0divT1 = []                                                                    # Ratio of inlet static temperatures [-]
    T1a = []                                                                        # Temperature just inside the blade [K]
    T1rmso = []                                                                     # Off design Root mean square of static temperature at inlet [K]
    P1rmso = []                                                                     # Off design root mean square of static pressure at the inlet [Pa]
    P1arms = []                                                                     # Total pressure just inside the bladed row [Pa]

    # Looping to find every variable of interest
    for i in range(0, len(Compressor.V0DivVcr)):
        eps.append(np.degrees(np.arctan((1 - BBF) * np.tan(np.radians(beta1rms[i])) / (1 + BBF * np.tan(np.radians(beta1rms[i])) ** 2))))     # (B40) Difference between compressor inlet relative flow angle and optimum incidence angle [deg]
        betaOpt.append(beta1rms[i] - eps[i])                                                # (B42) Optimum relative flow angle [deg]
        WL.append(W1rmso[i] * np.sin(np.radians(abs(betaOpt[i] - beta1rms[i]))))            # (B43) Component of relative velocity lost [m/s]
        dhInc.append((WL[i] ** 2) / 2)                                                      # (B44) Enthalpy loss due to incidence [J/kg]
        T1oRel.append(T1[i] + W1rmso[i] ** 2 / (2 * Fluid.Cp))                              # Off design inlet relative temperature [K]
        W1cr.append((2 * (Fluid.k - 1) / (Fluid.k + 1) * Fluid.R * T1oRel[i]) ** 0.5)       # Critical inlet relative velocity [m/s]    MSG: Do not understand this equation
        W1rmsEff.append(W1rmso[i] * np.cos(betaOpt[i] - beta1rms[i]))                       # Effective relative velocity [m/s]      MSG: Why is this the effective velocity?
        T0divT1.append(1 - (Fluid.k - 1) / (Fluid.k + 1) * (W1rmsEff[i] / W1cr[i]) ** 2)    # Ratio of inlet static to stagnation temperature [-]   MSG: Do not understand this equation
        T1a.append(T1oRel[i] * T0divT1[i])                                                  # Temperature just inside the blade [K]     MSG: Do not understand this equation
        T1rmso.append(T00i[i] - (Cm1rms[i] ** 2 / (2 * Fluid.Cp)))                          # Off design root mean square of static temperature at inlet [K], STAGNATION EQUATION 
        P1rmso.append(InletConditions.P00 * (T1rmso[i] / T00i[i]) ** (Fluid.k / (Fluid.k - 1)))     # Off design root mean square of static pressure at the inlet [Pa], ISENTROPIC PROPERTY
        P1arms.append(P1rmso[i] * np.exp((- dhInc[i]) / (T1a[i] * Fluid.R)))                # (B45) Total pressure just inside the bladed row [Pa]

    # Impeller work and losses
    U2o = U1to / (Compressor.r1 / Compressor.r2)                                            # Off design exit blade velocity [m/s], CONSTANT ANGULAR VELOCITY


    def density_impeller_outlet(rho2o):
        """ Calculate a series of velocities, temperatures and enthalpies corresponding to the impeller outlet density rho2o. 
            Finally, update density based on these properties. Function used for iteration """ 
        
        Vm2m = mdoto / (np.pi * rho2o * (2 * Compressor.r2) * b2)                               # (B49) Meridional component of exit absolute velocity [m/s]  , MASS BALANCE                
        VSL = U2o * (1 - sigma)                                                                 # (B51) Slip velocity [m/s]  , SLIP RELATION    MSG: Check that this chooses the correct slip factor      
        Vtheta2 = (U2o - Vm2m * np.tan(np.radians(beta2B)) - VSL)                               # (B50) Tangential component of exit absolute velocity [m/s]    , VELOCITY TRIANGLE    MSG: I think the sign before VM2m is correct here, but not sure
        T1orelrms = T1 + W1rmso ** 2 / (2 * Fluid.Cp)                                           # Relative root mean square temperature [K]   , STAGNATION TEMPERATURE
        
        T2orel = T1orelrms + ((U2o ** 2 - U1rmso ** 2) / (2 * Fluid.Cp))                        # (B52) Exit temperature in the relative reference frame [K]  , STAGNATION RELATION 
        Wtheta2 = U2o - Vtheta2                                                                 # (B53) Tangential component of relative exit velocity [m/s], VELOCITY TRIANGLE
        W2 = ((Vm2m ** 2) + (Wtheta2 ** 2)) ** 0.5                                              # (B54) Relative exit velocity [m/s]   ,  VELOCITY TRIANGLE
        T2o = (T2orel - ((W2 ** 2) / (2 * Fluid.Cp)))                                           # (B55) Off design point exit temperature [K]  , STAGNATION RELATION
        V2 = ((Vm2m ** 2) + (Vtheta2 ** 2)) ** 0.5                                              # (B56) Off design point absolute exit velocity [m/s]  ,  VELOCITY TRIANGLE  
        T2oabs = (T2o + (V2 ** 2) / (2 * Fluid.Cp))                                             # (B57) Off design point exit temperature in the absolute reference frame [K]   , STAGNATION RELATION
        dhaero = (Fluid.Cp * InletConditions.T00 * (T2oabs / InletConditions.T00 - 1))          # (B61) Aerodynamic enthalpy rise [J/kg]  , Cp*dT
        qaero = dhaero / (U2o ** 2)                                                             # (B60) Dimensionless actual head [-]   , WORK COEFFICIENT
        Df = (1 - W2 / W1to + (Compressor.kBL * qaero) / ((W1to / Compressor.U2) * ((ZB / np.pi) * (1 - 2 * Compressor.r1 / Compressor.D2) + 2 * 2 * Compressor.r1 / Compressor.D2)))   # (B59) Diffusion factor [-]  , EMPIRICAL
        dhBL = (0.05 * Df ** 2 * U2o ** 2)                                                      # (B58) Work loss due to blade loading [J/kg]   , ----
        Re = U2o * Compressor.D2 * rho1rms / Fluid.viscosity                                    # (B63) Reynolds number of the exit flow [-]   ,  ----
        dhDF = (0.01356 * rho2o * U2o ** 3 * Compressor.D2 ** 2 / (mdoto * Re ** 0.2))          # (B62) Impeller disk friction loss [J/kg]   , ----
        D1rms = np.sqrt((((2 * Compressor.r1) ** 2) + ((2 * Compressor.rh) ** 2)) / 2)          # Rootmean square of the diameter [m]  , RMS
        LenDivDia = 0.5 * (1 - (D1rms / 0.3048)) / (np.cos(np.radians(beta2B)))                 # (B65) Blade length to diameter ratio [-]   ,  ----
        #HydDiaDivExitDia = 1 / (ZB / (np.pi * np.cos(np.radians(beta2B)) + Compressor.D2 / b2)) + (2  * Compressor.r1 / Compressor.D2) / (2 / (1 - Fluid.k) + 2 * ZB / (np.pi * (1 + Fluid.k)) * np.sqrt(1 + (np.tan(np.radians(Compressor.beta1) ** 2) * (1 + Fluid.k ** 2 / 2))))     # (B66) Ratio of hydraulic diameter and exit diameter [-] MSG: One paranthesis wrong and use inducer hub-tip diameter ratio instead of specific heat ratio according to (B66)
        HydDiaDivExitDia = 1 / (ZB / (np.pi * np.cos(np.radians(beta2B))) + Compressor.D2 / b2) + (2  * Compressor.r1 / Compressor.D2) / (2 / (1 - Compressor.rh / Compressor.r1) + 2 * ZB / (np.pi * (1 + Compressor.rh / Compressor.r1)) * np.sqrt(1 + (np.tan(np.radians(Compressor.beta1) ** 2) * (1 + (Compressor.rh / Compressor.r1) ** 2 / 2))))     # (B66) Ratio of hydraulic diameter and exit diameter [-] MSG: This may be wrong, see additional_comments.tex
        WRelDivWext2 = 0.5 * ((Cm1rms / U2o) ** 2 + (D1rms / Compressor.D2) ** 2 + (W2 / W1to) ** 2 * ((Cm1rms / U2o) ** 2 + (2 * Compressor.r1 / Compressor.D2) ** 2))      # (B67) Ratio of mean relative velocity and impeller exit velocity^2 [-]
        dhSF = ((Compressor.kSF * Compressor.Cf * LenDivDia * WRelDivWext2 * U2o ** 2) / HydDiaDivExitDia)       # (B64) Skin friction loss [J/kg]
        #dhSF = Compressor.kSF * Compressor.Cf * LenDivDia / HydDiaDivExitDia * W1rmso ** 2     # (B64) MSG: Alternative

        dhid = (dhaero - dhInc - dhSF - dhDF - dhBL)                                            # (B68) Ideal enthalpy rise [J/kg]
        etaR = dhid / dhaero                                                                    # (B69) Impeller efficiency [-]  
        
        P2oabs = (P1arms * (etaR * dhaero / (Fluid.Cp * InletConditions.T00) + 1) ** (Fluid.k / (Fluid.k - 1))) # (B70) Iteration of the off design exit absolute pressure [Pa]
        
        # Check for invalid temperatures 
        for l in range(len(T2o)):
            if T2o[l] < 0:         
                raise ValueError('Temperature T2o is negative ' + str(T2o[l]) + '. Change range of V0/Vcr')
            
        # Find impeller outlet pressure and density
        P2o = (P2oabs / ((T2oabs / T2o) ** (Fluid.k / (Fluid.k - 1))))                          # (B71) Iteration of the off design exit pressure [Pa]      
        rho2oit = P2o / (Fluid.R * T2o)                                                         # (B72) Iteration of the off design exit density [kg/m^3]

        return [rho2oit, T2o, dhaero, dhBL, dhDF, dhSF, dhid, T2oabs, P2oabs, Vtheta2, Vm2m, Df, P2o, V2]


    def density_iteration():
        """ Determine impeller outlet density for varying mass flow by making an initial guess for the density and iterating it until it converges """

        dhEstim = U2o ** 2                                                                      # (B46) Initial approximation of enthalpy rise in impeller [J/kg], EULER EQUATION APPROXIMATION zero inlet swirl
        T2oEstAbs = (dhEstim / (Fluid.Cp * InletConditions.T00) + 1) * InletConditions.T00      # (B47) Estimate of the off design impeller exit total temperature [K], ENERGY/EULER EQUATION
        rho2o = Compressor.rho1 * (T2oEstAbs / InletConditions.T00) ** (1 / (Fluid.k - 1))      # (B48) Initial guess for impeller exit density [kg/m^3], ISENTROPIC RELATION
        rhoinit = rho2o    
        rho = [] 
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
        V2 = []

        # Determine impeller outlet density for different mass flows
        for i in range(0, len(Compressor.V0DivVcr)):
            #print('V0/Vcr =', Compressor.V0DivVcr[i])
            rhoafter = density_impeller_outlet(rho2o)[0][i]                
            if np.isnan(rhoafter):
                raise ValueError('Density calculation failed. Change range of V0/Vcr')
            rho_iterated = [rhoinit, rhoafter]                              

            # Iterate impeller outlet density until it converges
            while abs(((rho_iterated[- 1]) - (rho_iterated[- 2])) / rho_iterated[- 2]) > 0.0001:   
                rho2oiti, T2oi, dhaeroi, dhBLi, dhDFi, dhSFi, dhidi, T2oabsi, P2oabsi, Vtheta2i, Vm2mi, Dfi, P2oi, V2i = density_impeller_outlet(rho_iterated[- 1])
                rho_iterated.append(rho2oiti[i])                    
            
            rho_iterated.append(rho2oiti[i])    # MSG: If this fails, then V/Vcr is too small
            rho.append(rho2oiti[i])
            T2o.append(T2oi[i])
            dhaero.append(dhaeroi[i])
            dhBL.append(dhBLi[i])
            dhDF.append(dhDFi[i])
            dhSF.append(dhSFi[i])
            dhid.append(dhidi[i])
            T2oabs.append(T2oabsi[i])
            P2oabs.append(P2oabsi[i])
            Vtheta2.append(Vtheta2i[i])
            Vm2m.append(Vm2mi[i])
            Df.append(Dfi[i])
            P2o.append(P2oi[i])
            V2.append(V2i[i])
            
        return [rho, T2o, dhaero, dhBL, dhDF, dhSF, dhid, T2oabs, P2oabs, Vtheta2, Vm2m, Df, P2o, V2]

    rho2, T2o, dhaero, dhBL, dhDF, dhSF, dhid, T2oabs, P2oabs, Vtheta2, Vm2m, Df, P2o, V2 = density_iteration()

    rho2 = np.array(rho2)
    T2o = np.array(T2o)

    # Enthalpy loss and power contributions 
    dhaero = np.array(dhaero)                       # Aerodynamic enthalpy rise
    dhBL = np.array(dhBL)                           # Blade loading
    dhDF = np.array(dhDF)                           # Disk friction
    dhSF = np.array(dhSF)                           # Skin friciton
    dhid = np.array(dhid)                           # Ideal enthaply rise
    T2oabs = np.array(T2oabs)                       # Stagnation temperature at exit
    P2oabs = np.array(P2oabs)                       # Stagnation pressure at exit
    Vtheta2 = np.array(Vtheta2)                     # Tangential velocity component at exit
    Vm2m = np.array(Vm2m)                           # Meridonal velocity component at exit
    V2 = np.array(V2)                               # Absolute exit velocity
    Df = np.array(Df)                               # Diffusion factor
    P2o = np.array(P2o)                             # Stagnation pressure at exit after iteration
    C2o = []                                        # Absolute flow velocity at oulet

    # Finding impeller exit velocity
    for i in range(0, len(Compressor.V0DivVcr)):
        C2o.append((Vm2m[i] ** 2 + Vtheta2[i] ** 2) ** 0.5)
    M2o = C2o / np.sqrt(Fluid.k * Fluid.R * T2oabs)         # MSG: Think this is wrong. Here, speed of sound calculation uses absolute/stagnation temperature, but should use static temperature. Does not affect other properties

    # Recirculation loss, work done on fluid going back 
    dhRC = []
    alpha2 = []
    for i in range(0, len(Compressor.V0DivVcr)):
        alpha2.append(np.degrees(np.arctan(Vtheta2[i] / Vm2m[i])) )                         # Exit velocity flow angle [deg]
        dhRC.append(0.02 * np.sqrt(np.tan(np.radians(alpha2[i]))) * Df[i] ** 2 * U2o ** 2)  # (B73) Enthalpy loss from recirculation [J/kg]

    # Diffuser calculations and exit losses
    etad = Compressor.etad                          # Diffuser efficiency
    CpDi = 1 - (Compressor.AR ** (- 2))             # Ideal pressure recovery coefficient 
    CpD = etad * CpDi                               # Pressure recovery coefficient
    M3o = []                                        # Mach number
    P3oabs = []                                     
    dhVLD = []                                      # Vaneless diffuser loss
    P3o = []                                        # Outlet stagnation pressure
    C3o = []                                        # Outlet velocity
    P03o = []                                       
    for i in range(0, len(Compressor.V0DivVcr)):
        P3o.append(P2o[i] + CpD * 0.5 * rho2[i] * C2o[i] ** 2)
        C3o.append(C2o[i] / Compressor.AR)
        #T3o = T3oabs - C3o ** 2 / (2 * Fluid.Cp)
        
        M3o.append(C3o[i] / (np.sqrt(Fluid.k * Fluid.R * T2oabs[i])))       # MSG: Think this is wrong. Here, speed of sound calculation uses absolute/stagnation temperature, but should use static temperature. It is also used at the wrong point. Should be T3o. Any way to obtain T3o/T3oabs?
        P3oabs.append(P3o[i] * (1 + (Fluid.k - 1) / 2 * M3o[i] ** 2) ** (Fluid.k / (Fluid.k - 1)))    # (B85) Absolute diffuser throat pressure [Pa]      
        dhVLD.append(Fluid.Cp * T2oabs[i] * ((P3o[i] / P3oabs[i]) ** ((Fluid.k - 1) / Fluid.k) - (P3o[i] / P2oabs[i]) ** ((Fluid.k - 1) / Fluid.k)))  # (B86) Vaneless diffuser loss [J/kg]
        
    # Overall performance
    etao = []
    etaoAlt = []
    Pro = []
    
    for i in range(0, len(Compressor.V0DivVcr)):
        etaAppend = (dhaero[i] - (dhInc[i] + dhBL[i] + dhSF[i] + dhVLD[i])) / (dhaero[i] + dhRC[i] + dhDF[i])       # MSG: Why put dhRC and dhDF in the denominator and not in nominator?
        etao.append(etaAppend)      
        P3append = P3o[i] 
        # P3append = P3oabs[i]      
        Pro.append(P3append / InletConditions.P00)
    
    enthalpy_rise_loss = [dhaero, dhInc, dhBL, dhSF, dhVLD, dhRC, dhDF]
    return Pro, P3o, T2oabs, mdoto, etao, U2o, M2o, enthalpy_rise_loss, V2, etaoAlt


def off_design_performance(Compressor, Fluid, InletConditions, IterationMatrix):
    """ Calculate off-design performance for varying rotational speed by calling off_design_performance_rpm function """

    results_off_design = []
    N_errors = 0
    for iN in range(0, len(Compressor.N_off_design_arr)):
        print('\nRotational velocity =', Compressor.N_off_design_arr[iN], 'RPM')
        Compressor.V0DivVcr = np.linspace(Compressor.V0DivVcr_off_design_arr[2 * iN], Compressor.V0DivVcr_off_design_arr[2 * iN + 1], 400)
        Pro, P3o, T2oabs, mdoto, etao, U2o, M2o, enthalpy_rise_loss, V2, etaoAlt = off_design_performance_rpm(Compressor.N_off_design_arr[iN], Compressor, Fluid, InletConditions, IterationMatrix)

        # Check if all efficiencies are in the range [0, 1]
        etao = np.array(etao)
        if any(etao < 0) or any(etao > 1):  
            print(f'\tError: Rotational speed = {Compressor.N_off_design_arr[iN]} RPM has efficiency outside range [0, 1]. Efficiencies =', etao)
            N_errors += 1
        
        # Make arrays for corrected mass flow and corrected speed using reference conditions
        T_ref = 288.15
        P_ref = 1.01300e5 
        correctedMass = mdoto * (InletConditions.P00 / P_ref) / np.sqrt(InletConditions.T00 / T_ref)       
        correctedSpeed = Compressor.N_off_design_arr[iN] / np.sqrt(InletConditions.T00 / T_ref)

        results_off_design.append({
            'N': Compressor.N_off_design_arr[iN],
            'Pro': np.array(Pro),
            'P3o': np.array(P3o),
            'T2oabs': np.array(T2oabs),
            'mdoto': np.array(mdoto),
            'etao': np.array(etao),
            'U2o': np.array(U2o),
            'M2o': np.array(M2o),
            'mdoto_corrected': np.array(correctedMass),
            'N_corrected': np.array(correctedSpeed),
            'enthalpy_rise_loss': np.array(enthalpy_rise_loss),
            'V2': np.array(V2),
            'etaoAlt': np.array(etaoAlt)
        })

    print('\nOff-design_performance calculated successfully.')
    print(f'Number of errors in efficiency: {N_errors} out of {len(Compressor.N_off_design_arr)} rotational speeds')
        
    return results_off_design


if __name__ == '__main__':
    pass
