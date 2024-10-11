"""
The following Python code uses the inducer and impeller geometry from the geometry code
(geometry.py) and predicts the performance of the compressor at off-design points.

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


References: Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244
            Galvas, Michael R. (1973). Fortran program for predicting off-design performance of centrifugal compressors, https://ntrs.nasa.gov/citations/19740001912 

Author(s): Petter Resell (summer intern, 2024), Martin Spillum Grønli (SINTEF Energy Research, 2024)
"""


### Import---------------------------------------------------------------------------------------------------       
import numpy as np
import matplotlib.pyplot as plt


def off_design_performance_rpm(Nrpm, Compressor, Fluid, InletConditions, IterationMatrix):
    """ Function takes only rpm and decides all flow properties. It is mainly based on iteration of density. """
    """ ---------- Control, Determining point of interest ---------- """
    bladeAngle = Compressor.bladeAngle                                   # Blade angle of interest MSG: This must be updated after geometry calculations. Need to pick one
    bladeNumber = Compressor.bladeNumber                                 # Blade number of interest MSG: This must be updated after geometry calculations. Need to pick one

    iB = int(np.where(IterationMatrix.beta2BArr == bladeAngle)[0])                   # Index for blade angle
    iZ = int(np.where(IterationMatrix.ZBarr == bladeNumber)[0])                      # Index for blade number
    beta2B = IterationMatrix.beta2BArr[iB]                                           # 
    ZB = IterationMatrix.ZBarr[iZ]                                                   # 
    b2 = IterationMatrix.b2Mat[iB][iZ]                                               # Compressor  outlet height 
    sigma = IterationMatrix.sigmaMat[iB][iZ]                                         # Slip factor        

    if np.isnan(b2):
        raise Exception(f"Blade number and blade angle is not found to be among valid cases! \n Chose another point or adjust parameters!")

    """ ---------- Parameters, non-variable ---------- """
    VelCr = np.sqrt(2 * Fluid.k / (Fluid.k + 1) * Fluid.R * InletConditions.T00)    # Critical Velocity wrt resonance [m/s]  MSG: This has unit m/s sqrt(kg/mol). Also do not understand this equation.
    Cm1h = Compressor.V0DivVcr * VelCr                                              # Absolute meridional velocity of the hub [m/s]
    T00i = np.full(len(Compressor.V0DivVcr), InletConditions.T00)                   # Inlet temperature array [K]


    # ---------------------------- Calculation of Inlet Velocity Triangles and Compressor Weight Flow (Swirl free) ----------------------------
    curve1rms = np.sqrt((Compressor.curvet1 ** 2 + Compressor.curveh1 ** 2) / 2)    # Root mean square of the inducer inlet hub and tip wall curvature [m^-1]
    r1rms = np.sqrt((Compressor.r1 ** 2 + Compressor.rh1 ** 2) / 2)                 # Root mean square of the hub and tip inlet radius [m]  MSG: Why is RMS used here?

    # Geometrical inlet properties
    h0 = r1rms - Compressor.rh1                 # (B4) Spacing for numerical integration [-]
    h1 = Compressor.r1 - r1rms                  # (B5) Spacing for numerical integration [-]

    # Absolute meridonal velocities
    Cm1rms = Cm1h * np.exp((h0 / 2) * (Compressor.curveh1 + curve1rms))     # (B2) Absolute meridional root mean square velocity [m/s]
    Cm1t = Cm1h * np.exp(((h0 + h1) / 6) * ((2 - (h1 / h0)) * Compressor.curveh1 + (((h0 + h1) ** 2) / (h1 * h0)) * curve1rms + (2 - (h0 / h1)) * Compressor.curvet1))      # (B3) Absolute meridional velocity of the tip [m/s]

    # Velocities through the normal flow area at inlet
    Cm1hn = Cm1h * np.cos(np.radians(Compressor.x))                         # (B6) Normal component of the absolute hub velocity [m/s], through flow area
    Cm1tn = Cm1t * np.cos(np.radians(Compressor.x))                         # (B6) Normal component of the absolute tip velocity [m/s], through flow area
    Cm1rmsn = Cm1rms * np.cos(np.radians(Compressor.x))                     # (B6) Normal component of the root mean square velocity [m/s], through flow area

    # Densities from stagnation temperature and different velocites
    rho1h = InletConditions.rho0 * (1 - (Cm1h ** 2 / (2 * Fluid.Cp * InletConditions.T00))) ** (1 / (Fluid.k - 1))                # Inlet density of the fluid at the hub [kg/m^3], from stagnation temperature
    rho1rms = InletConditions.rho0 * (1 - (Cm1rms ** 2 / (2 * Fluid.Cp * InletConditions.T00))) ** (1 / (Fluid.k - 1))            # Root mean square density of the fluid [kg/m^2], 
    rho1t = InletConditions.rho0 * (1 - (Cm1t ** 2 / (2 * Fluid.Cp * InletConditions.T00))) ** (1 / (Fluid.k - 1))                # Inlet density of the fluid at the blade tip [kg/m^3]

    mdoto = 2 * np.pi * (((h0 + h1) / 6) * ((2 - (h1 / h0)) * (rho1h * Compressor.rh1 * Cm1hn) + ((h0 + h1) ** 2 / (h0 * h1)) * (rho1rms * r1rms * Cm1rmsn) + (2 - (h0 / h1)) * (rho1t * Compressor.r1 * Cm1tn)))      # (B7) Off Design point mass flow rate [kg/s]

    # Blade speeds
    U1to = np.pi * Compressor.No * Nrpm * 2 * Compressor.r1 / 60            # Off design blade tip velocity [m/s]
    U1ho = np.pi * Compressor.No * Nrpm * 2 * Compressor.rh1 / 60           # Off design blade hub velocity [m/s]
    U1rmso = ((U1to ** 2 + U1ho ** 2) / 2) ** 0.5                           # Off design root mean square velocity [m/s]

    # Relative velocities at inlet
    W1ho = (Cm1h ** 2 + U1ho ** 2) ** 0.5                           # Off design hub relative velocity [m/s]  
    W1to = (Cm1t ** 2 + U1to ** 2) ** 0.5                           # Off design tip relative velocity [m/s]  
    W1rmso = ((W1ho ** 2 + W1to ** 2) / 2) ** 0.5                   # Root mean square relative velocity [m/s]

    # Relative flow inlet angles beta (not blade angle, but flow angle)
    beta1t = []                                 # Hub inlet relative angle [degrees]
    beta1h = []                                 # Tip inlet relative angle [degrees]
    beta1rms = []                               # Root mean square inlet relative angle [degrees]

    T1 = T00i - (Cm1rms ** 2) / (2 * Fluid.Cp)  # Inlet static temperature [K]
    for i in range (0, len(Compressor.V0DivVcr)):
        #beta1t.append(np.degrees(np.arctan(Cm1t[i] / U1to)))               # Tip inlet relative angle [degrees]    MSG: Should this be U1to/Cm1t instead? Seems to be different angle than what is used when finding eps. See https://ntrs.nasa.gov/api/citations/19930083715/downloads/19930083715.pdf. Also, is the tangential velocity of the fluid at the impeller tip inlet upstream of the blade assumed to be zero? Ref. eq. B26 and B25
        beta1t.append(np.degrees(np.arctan(U1to / Cm1t[i])))
        #beta1h.append(np.degrees(np.arctan(Cm1h[i] / U1ho)))               # Hub inlet relative angle [degrees]    MSG: Should this be U1to/Cm1t instead? Seems to be different angle than what is used when finding eps. See https://ntrs.nasa.gov/api/citations/19930083715/downloads/19930083715.pdf. Also, is the tangential velocity of the fluid at the impeller tip inlet upstream of the blade assumed to be zero? Ref. eq. B26 and B25
        beta1h.append(np.degrees(np.arctan(U1ho / Cm1h[i])))

        beta1rms.append(((beta1t[i] ** 2 + beta1h[i] ** 2) / 2) ** 0.5)     # Root mean square inlet relative angle [degrees]



    ###---------------------- Inducer Incidence Loss ----------------------

    BBF = 1 - (ZB * Compressor.tu) / (2 * np.pi * r1rms)     # (41) Blade Blockage Factor [-]        
    # Declaring variables before use
    eps = []                                        # Difference between compressor inlet relative flow angle and optimum incidence angle [degrees]
    betaOpt = []                                    # Optimum relative flow angle [degrees]
    WL = []                                         # Relative velocity loss [m/s]
    dhInc = []                                      # Enthalpy loss due to incidence [J/kg]
    T1oRel = []                                     # Off design inlet relative temperature [K]
    W1cr = []                                       # Critical inlet relative velocity [m/s]
    W1rmsEff = []                                   # Effective relative velocity [m/s]
    T0divT1 = []                                    # Ratio of inlet static temperatures [-]
    T1a = []                                        # Temperature just inside the blade [K]
    T1rmso = []                                     # Off design Root mean square of static temperature at inlet [K]
    P1rmso = []                                     # Off design root mean square of static pressure at the inlet [Pa]
    P1arms = []                                     # Total pressure just inside the bladed row [Pa]

    # Looping to find every variable of interest
    for i in range(0, len(Compressor.V0DivVcr)):
        eps.append(np.degrees(np.arctan((1 - BBF) * np.tan(np.radians(beta1rms[i])) / (1 + BBF * np.tan(np.radians(beta1rms[i])) ** 2))))     # (B40) Difference between compressor inlet relative flow angle and optimum incidence angle [degrees]
        betaOpt.append(beta1rms[i] - eps[i])                                        # (B42) Optimum relative flow angle [degrees]
        WL.append(W1rmso[i] * np.sin(np.radians(abs(betaOpt[i] - beta1rms[i]))))    # (B43) Component of relative velocity lost [m/s]
        
        dhInc.append((WL[i] ** 2) / 2)                                              # (B44) Enthalpy loss due to incidence [J/kg]
        
        T1oRel.append(T1[i] + W1rmso[i] ** 2 / (2 * Fluid.Cp))                              # Off design inlet relative temperature [K]
        W1cr.append((2 * (Fluid.k - 1) / (Fluid.k + 1) * Fluid.R * T1oRel[i]) ** 0.5)       # Critical inlet relative velocity [m/s]    MSG: Do not understand this equation
        W1rmsEff.append(W1rmso[i] * np.cos(betaOpt[i] - beta1rms[i]))                       # Effective relative velocity [m/s]      MSG: Why is this the effective velocity?
        T0divT1.append(1 - (Fluid.k - 1) / (Fluid.k + 1) * (W1rmsEff[i] / W1cr[i]) ** 2)    # Ratio of inlet static to stagnation temperature [-]   MSG: Do not understand this equation
        T1a.append(T1oRel[i] * T0divT1[i])                                                  # Temperature just inside the blade [K]     MSG: Do not understand this equation
        T1rmso.append(T00i[i] - (Cm1rms[i] ** 2 / (2 * Fluid.Cp)))                          # Off design root mean square of static temperature at inlet [K], STAGNATION EQUATION 
        P1rmso.append(InletConditions.P00 * (T1rmso[i] / T00i[i]) ** (Fluid.k / (Fluid.k - 1)))     # Off design root mean square of static pressure at the inlet [Pa], ISENTROPIC PROPERTY
        P1arms.append(P1rmso[i] * np.exp((- 1 * dhInc[i]) / (T1a[i] * Fluid.R)))            # (B45) Total pressure just inside the bladed row [Pa]



    ### ---------------------- Impeller Work and Losses ----------------------
    U2o = U1to / (Compressor.r1 / Compressor.r2)                # Off design exit blade velocity [m/s], CONSTANT ANGULAR VELOCITY
    dhEstim = U2o ** 2                                          # (B46) Initial approximation of enthalpy rise in impeller [J/kg], EULER EQUATION APPROXIMATION zero inlet swirl
    T2oEstAbs = (dhEstim / (Fluid.Cp * InletConditions.T00) + 1) * InletConditions.T00      # (B47) Estimate of the off design impeller exit total temperature [K], ENERGY/EULER EQUATION
    rho2o = Compressor.rho1 * (T2oEstAbs / InletConditions.T00) ** (1 / (Fluid.k - 1))      # (B48) Off design impeller exit density [kg/m^3], ISENTROPIC RELATION


    def Densityiteration(rho2o):
        """
        The Densityiteration(rho2o) function takes an initial guess of the impeller outlet density using
        Equation (B48) above. It then uses this initial guess to calculate a series of velocities,
        temperatures and enthalpies corresponding to this initial guess before re-calculating the
        density.
        """ 
        Vm2m = mdoto / (np.pi * rho2o * (2 * Compressor.r2) * b2)                       # (B49) Meridional component of exit absolute velocity [m/s]  , MASS BALANCE                
        VSL = U2o * (1 - sigma)                                                         # (B51) Slip velocity [m/s]  , SLIP RELATION    MSG: Check that this chooses the correct slip factor      
        Vtheta2 = (U2o - Vm2m * np.tan(np.radians(- beta2B)) - VSL)                     # (B50) Tangential component of exit absolute velocity [m/s]    , VELOCITY TRIANGLE    MSG: I think the sign before VM2m is correct here, but not sure
        T1orelrms = T1 + W1rmso ** 2 / (2 * Fluid.Cp)                                   # Relative root mean square temperature [K]   , STAGNATION TEMPERATURE
        
        T2orel = T1orelrms + ((U2o ** 2 - U1rmso ** 2) / (2 * Fluid.Cp))                # (B52) Exit temperature in the relative reference frame [K]  , STAGNATION RELATION    MSG: Something wrong here. Here and below we have multiple equations for the same thing. This line I wrote myself, and I think this is correct
        #T2orel = T1orelrms + ((U2o ** 2 - Compressor.U1t ** 2) / (2 * Fluid.Cp))       # (B52) Exit temperature in the relative reference frame [K]  , STAGNATION RELATION    MSG: Something wrong here. Here and below we have multiple equations for the same thing.
        #T2orel = T1orelrms + ((U2o ** 2 - U1to ** 2) / (2 * Fluid.Cp))                 # (B52) Exit temperature in the relative reference frame [K]  , STAGNATION RELATION    MSG: Something wrong here. Here and below we have multiple equations for the same thing.
        #T2orel = T1 + ((U2o ** 2 - U1t ** 2) / (2 * Cp))                               # (B52) Exit temperature in the relative reference frame [K]  , STAGNATION RELATION    MSG: Something wrong here. Here and below we have multiple equations for the same thing.
        
        Wtheta2 = U2o - Vtheta2                                                         # (B53) Tangential component of relative exit velocity [m/s], VELOCITY TRIANGLE
        W2 = ((Vm2m ** 2) + (Wtheta2 ** 2)) ** 0.5                                      # (B54) Relative exit velocity [m/s]   ,  VELOCITY TRIANGLE
        T2o = (T2orel - ((W2 ** 2) / (2 * Fluid.Cp)))                                   # (B55) Off design point exit temperature [K]  , STAGNATION RELATION
        V2 = ((Vm2m ** 2) + (Vtheta2 ** 2)) ** 0.5                                      # (B56) Off design point absolute exit velocity [m/s]  ,  VELOCITY TRIANGLE  
        T2oabs = (T2o + (V2 ** 2) / (2 * Fluid.Cp))                                     # (B57) Off design point exit temperature in the absolute reference frame [K]   , STAGNATION RELATION
        dhaero = (Fluid.Cp * InletConditions.T00 * (T2oabs / InletConditions.T00 - 1))  # (B61) Aerodynamic enthalpy rise [J/kg]  , Cp*dT
        qaero = dhaero / (U2o ** 2)                                                     # (B60) Dimensionless actual head [-]   , WORK COEFFICIENT
        Df = (1 - W2 / W1to + (Compressor.kBL * qaero) / ((W1to / Compressor.U2) * ((ZB / np.pi) * (1 - 2 * Compressor.r1 / Compressor.D2) + 2 * 2 * Compressor.r1 / Compressor.D2)))   # (B59) Diffusion factor [-]  , EMPIRICAL
        dhBL = (0.05 * Df ** 2 * U2o ** 2)                                                      # (B58) Work loss due to blade loading [J/kg]   , ----
        Re = U2o * Compressor.D2 * rho1rms / Fluid.viscosity                                    # (B63) Reynolds number of the exit flow [-]   ,  ----
        dhDF = (0.01356 * rho2o * U2o ** 3 * Compressor.D2 ** 2 / (mdoto * Re ** 0.2))          # (B62) Impeller disk friction loss [J/kg]   , ----
        D1rms = np.sqrt((((2 * Compressor.r1) ** 2) + ((2 * Compressor.rh1) ** 2)) / 2)         # Rootmean square of the diameter [m]  , RMS
        LenDivDia = 0.5 * (1 - (D1rms / 0.3048)) / (np.cos(np.radians(beta2B)))                 # (B65) Blade length to diameter ratio [-]   ,  ----
        HydDiaDivExitDia = 1 / (ZB / (np.pi * np.cos(np.radians(beta2B)) + Compressor.D2 / b2)) + (2  * Compressor.r1 / Compressor.D2) / (2 / (1 - Fluid.k) + 2 * ZB / (np.pi * (1 + Fluid.k)) * np.sqrt(1 + (np.tan(np.radians(Compressor.beta1) ** 2) * (1 + Fluid.k ** 2 / 2))))     # (B66) Ratio of hydraulic diameter and exit diameter [-]
        WRelDivWext = 0.5 * ((Cm1rms / U2o) ** 2 + (D1rms / Compressor.D2) ** 2 + (W2 / W1to) ** 2 * ((Cm1rms / U2o) ** 2 + (2 * Compressor.r1 / Compressor.D2) ** 2))      # (B67) Ratio of mean relative velocity and impeller exit velocity^2 [-]
        dhSF = ((Compressor.kSF * Compressor.Cf * LenDivDia * WRelDivWext * U2o ** 2) / HydDiaDivExitDia)       # (B64) Skin friction loss [J/kg]
        dhid = (dhaero - dhInc - dhSF - dhDF - dhBL)                                                            # (B68) Ideal enthalpy rise [J/kg]
        etaR = dhid / dhaero                                                                                    # (B69) Impeller efficiency [-]
        P2oabs = (P1arms * (etaR * dhaero / (Fluid.Cp * InletConditions.T00) + 1) ** (Fluid.k / (Fluid.k - 1))) # (B70) Iteration of the off design exit absolute pressure [Pa]
        """ Checking for invalid temperatures """
        for l in range(len(T2o)):
            if T2o[l] < 0:             # MSG: Added this to raise error if negative temperatures
                #plt.show()
                raise ValueError('Temperature T2o is negative ' + str(T2o[l]) + '. Change range of V0DivVcr in settings.py') 
            
        """ Finding impeller outlet temp. and density"""
        P2o = (P2oabs / ((T2oabs / T2o) ** (Fluid.k / (Fluid.k - 1))))                  # (B71) Iteration of the off design exit pressure [Pa]      
        rho2oit = P2o / (Fluid.R * T2o)                                                 # (B72) Iteration of the off design exit density [kg/m^3]

        return [rho2oit, T2o, dhaero, dhBL, dhDF, dhSF, dhid, T2oabs, P2oabs, Vtheta2, Vm2m, Df, P2o]


    def Density():
        """
        The Density() function selects the value corresponding to the variables listed below when the
        output density from the DensityIteration(rho2o) function is within 0.1% of the input. I.e. iterating until the output density stabilizes.
        """
        rhoinit = rho2o
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

        """ Iterating through the length of array V0DivVcr that is the ratio of the inlet velocity to the critical velocity """
        for i in range(0, len(Compressor.V0DivVcr)):
            print('\n\n', i)
            rhoafter = Densityiteration(rho2o)[0][i]                # Initial density. MSG: Does i represent a given mass flow? Inlet velocity (set in settings.py by self.V0DivVcr) controls the mass flow rate
            rho = [rhoinit, rhoafter]                               # Control array for iteration, reset for each index i

            while abs(((rho[- 1]) - (rho[- 2])) / rho[- 2]) > 0.0001:
                rho2oiti, T2oi, dhaeroi, dhBLi, dhDFi, dhSFi, dhidi, T2oabsi, P2oabsi, Vtheta2i, Vm2mi, Dfi, P2oi = Densityiteration(rho[- 1])
                print('rho2oiti', rho2oiti)
                rho.append(rho2oiti[i])
                    
            rho.append(rho2oiti[i])
            RHO.append(rho2oiti[i])
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
            
        return [RHO, T2o, dhaero, dhBL, dhDF, dhSF, dhid, T2oabs, P2oabs, Vtheta2, Vm2m, Df, P2o]
        
        # MSG: Old version of iteration below where I do not understand the purpose of having two conditions (0.000001 and 0.001)
        """
        #Iterating through the length of array V0DivVcr that is the ratio of the inlet velocity to the critical velocity 
        for i in range(0, len(Compressor.V0DivVcr)):
            rhoafter = Densityiteration(rho2o)[0][i]                # Initial density. MSG: Does i represent a given mass flow? Inlet velocity controls the mass flow rate
            rho = [rhoinit, rhoafter]                               # Control array for iteration, reset for each index i

            while abs((rho[- 1]) - (rho[- 2])) > 0.0000001:
                print(i)
                # If the difference between densities at the current and previous iteration is greater than one thousandth the density is added to the array but iteration continues. 
                #        If the difference between iterated densities is less than one thousandth then the density at index i is replaced and all other 
                #            properties at the same index is found

                rho2oiti, T2oi, dhaeroi, dhBLi, dhDFi, dhSFi, dhidi, T2oabsi, P2oabsi, Vtheta2i, Vm2mi, Dfi, P2oi = Densityiteration(rho[- 1])
                if abs((rho[- 1]) - (rho[- 2])) > 0.001:
                    print('\n Third', len(RHO), i)
                    rho.append(rho2oiti[i])
                    
                else:
                    
                    if len(RHO) <= i:
                        print('\n First', len(RHO), i)
                        rho.append(rho2oiti[i])
                        RHO.append(rho2oiti[i])
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
                    else:
                        print('\n Second', len(RHO), i)
                        rho.append(rho2oiti[i])
                        
        return [RHO, T2o, dhaero, dhBL, dhDF, dhSF, dhid, T2oabs, P2oabs, Vtheta2, Vm2m, Df, P2o, ]
        """

    rho2, T2o, dhaero, dhBL, dhDF, dhSF, dhid, T2oabs, P2oabs, Vtheta2, Vm2m, Df, P2o = Density()

    
    rho2 = np.array(rho2)
    T2o = np.array(T2o)

    # Enthalpy loss and power contributions 
    dhaero = np.array(dhaero)                   # Aerodynamic enthalpy rise
    dhBL = np.array(dhBL)                       # Blade loading
    dhDF = np.array(dhDF)                       # Disk friction
    dhSF = np.array(dhSF)                       # Skin friciton
    dhid = np.array(dhid)                       # Ideal enthaply rise
    T2oabs = np.array(T2oabs)                   # Stagnation temperature at exit
    P2oabs = np.array(P2oabs)                   # Stagnation pressure at exit
    Vtheta2 = np.array(Vtheta2)                 # Tangential velocity component at exit
    Vm2m = np.array(Vm2m)                       # Meridonal velocity component at exit
    Df = np.array(Df)                           # Diffusion factor
    P2o = np.array(P2o)                         # Stagnation pressure at exit after iteration
    C2o = []                                    # Absolute flow velocity at oulet
  

    """ Finding impeller exit velocity """
    for i in range(0, len(Compressor.V0DivVcr)):
        #print(Vm2m, Vtheta2)
        C2o.append((Vm2m[i] ** 2 + Vtheta2[i] ** 2) ** 0.5)
    M2o = C2o / np.sqrt(Fluid.k * Fluid.R * T2oabs)


    ### Recirculation loss, work done on fluid going back 
    dhRC = []
    alpha2 = []
    for i in range(0, len(Compressor.V0DivVcr)):
        alpha2.append(np.degrees(np.arctan(Vtheta2[i] / Vm2m[i])) )                         # Exit velocity flow angle [degrees]
        dhRC.append(0.02 * np.sqrt(np.tan(np.radians(alpha2[i]))) * Df[i] ** 2 * U2o ** 2)  # (B73) Enthalpy loss from recirculation [kJ/kg]


    ### Exit losses--------------------------------------------------------------------------------------------
    """ Diffuser calculations """
    etad = Compressor.etad                          # (...) TODO
    CpDi = 1 - (Compressor.AR ** (- 2))             # (...) TODO
    CpD = etad * CpDi                               # (...) TODO
    M3o = []                                        # Mach number
    P3oabs = []                                     # Stagnation outlet pressure
    dhVLD = []                                      # vane less diffuser loss
    P3o = []                                        # Outlet static pressure
    C3o = []                                        # Outlet velocity
    P03o = []                                       # (...) TODO
    for i in range(0, len(Compressor.V0DivVcr)):
        P3o.append(CpD * 0.5 * rho2[i] * C2o[i] ** 2 + P2o[i])
        C3o.append(C2o[i] / Compressor.AR)
        P03o.append(P3o[i] + 0.5 * rho2[i] * C3o[i] ** 2)
        M3o.append(C3o[i] / (np.sqrt(Fluid.k * Fluid.R * T2oabs[i])))
        P3oabs.append(P3o[i] * (1 + (Fluid.k - 1) / 2 * M3o[i] ** 2) ** (Fluid.k / (Fluid.k - 1)))    # (B85) Absolute diffuser throat pressure [Pa]      
        dhVLD.append(Fluid.Cp * T2oabs[i] * ((P3o[i] / P3oabs[i]) ** ((Fluid.k - 1) / Fluid.k) - (P3o[i] / P2oabs[i]) ** ((Fluid.k - 1) / Fluid.k)))  # (B86) Vaneless diffuser loss [kJ/kg]
    

    ### Overall performance-----------------------------------------------------------------------------------
    etao = []
    Pro = []
    for i in range(0, len(Compressor.V0DivVcr)):
        etaAppend = (dhaero[i] - (dhInc[i] + dhBL[i] + dhSF[i] + dhVLD[i])) / (dhaero[i] + dhRC[i] + dhDF[i])    
        etao.append(etaAppend)      
        P3append = P3o[i]           
        # P3append = P3oabs[i]      
        Pro.append(P3append / InletConditions.P00)

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
    # plt.plot(mdoto, Pro, 'r', label = str(Nrpm) +" RPM", linewidth = 2)
    # plt.xlabel("Mass flow rate (kg/s)")
    # plt.ylabel("Pressure ratio (-)")
    # plt.title("Compressor Map")


    ### Matching-------------------------------------------------------------------------------------------
    # comp_power = [100]      # Compressor power [kW]
    # comp_power = [round(IterationMatrix.WxMat[iB][iZ], 2) ]      # Compressor power [kW]
    # # speed = [120]           # Shaft Speed [krpm]    MSG: Not in use
    # T1 = InletConditions.T00                # Inlet stagnation temperature [K]
    # cp = Fluid.Cp              # Specific heat at constant pressure [kJ]
    # y = Fluid.k                 # Ratio of specific heats
    # # P = 0.1013              # Inlet pressure [MPa]
    # eta = IterationMatrix.etaMat[iB][iZ]              # Assumed efficiency

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

    print('Current RPM: ' + str(round(Nrpm, 2)) + '     Inlet blade tip speed: ' + str(round(U1to, 2)) + ' m/s     Outlet blade tip speed: ' + str(round(U2o, 2)) + ' m/s' )
  
    return Pro, P03o, T2oabs, mdoto, etao, U2o, M2o


def off_design_performance(Compressor, Fluid, InletConditions, IterationMatrix):
    """ Iterate through different rotational speeds and calculate the off-design performance of the compressor for each speed by calling off_design_performance_rpm function. Plot the results. """

    # Running off-desing_performance.py for rotational speeds of 15000rpm to the critical speed with increments of 5000
    Narr = np.arange(30000, Compressor.Ncrit, 5000)    
    #Narr = np.arange(80000, 120000, 5000)    
    N_errors = 0
    
    constant_eff_lines = [0.82, 0.83, 0.84, 0.85, 0.86, 0.87]    # Constant efficiency lines
    colors = ['r', 'g', 'b', 'c', 'm', 'y']                      # Colors for the constant efficiency lines
    constant_eff_line_mdot = [[[[], []] for i in range(len(constant_eff_lines))] for _ in range(len(Narr))]    # Constant efficiency line 1
    constant_eff_line_pr = [[[[], []] for i in range(len(constant_eff_lines))] for _ in range(len(Narr))]      # Constant efficiency line 2
    

    # ------------------- Plotting from off_design_performance.py -------------------
    for iN in range(0, len(Narr)):                                                                      # Running and plotting for all rpm N's
        Pro, P03o, T2oabs, mdoto, etao, U2o, M2o = off_design_performance_rpm(Narr[iN], Compressor, Fluid, InletConditions, IterationMatrix)    # Running everything
        Pro = np.array(Pro)                                                                             # Pressure ratio
        P03o = np.array(P03o)                                                                           # Outlet pressure
        T2oabs = np.array(T2oabs)                                                                       # Outlet temperature
        mdoto = np.array(mdoto)                                                                         # Mass flow
        etao = np.array(etao)                                                                           # Efficiency
        correctedMass = mdoto * (P03o / InletConditions.P00) / np.sqrt(T2oabs / InletConditions.T00)    # Corrected mass flow parameter

        # Finding the constant efficiency lines. MSG: Need to connect some of the lines
        for i in range(len(constant_eff_lines)):
            if any(etao >= constant_eff_lines[i]):
                if etao[0] <= constant_eff_lines[i]:
                    index_0 = np.argmax(etao >= constant_eff_lines[i])
                    constant_eff_line_mdot[iN][i][0] = mdoto[index_0]
                    constant_eff_line_pr[iN][i][0] = Pro[index_0]
                    if etao[- 1] <= constant_eff_lines[i]:
                        index_1 = np.argmax(etao[index_0:] < constant_eff_lines[i])
                        constant_eff_line_mdot[iN][i][1] = mdoto[index_0 + index_1]
                        constant_eff_line_pr[iN][i][1] = Pro[index_0 + index_1]
                    else:
                        constant_eff_line_mdot[iN][i][1] = np.nan
                        constant_eff_line_pr[iN][i][1] = np.nan
                else:
                    if etao[- 1] <= constant_eff_lines[i]:
                        index_1 = np.argmax(etao < constant_eff_lines[i])
                        constant_eff_line_mdot[iN][i][1] = mdoto[index_1]
                        constant_eff_line_pr[iN][i][1] = Pro[index_1]
                        constant_eff_line_mdot[iN][i][0] = np.nan
                        constant_eff_line_pr[iN][i][0] = np.nan
                    else:
                        constant_eff_line_mdot[iN][i][0] = np.nan
                        constant_eff_line_pr[iN][i][0] = np.nan
                        constant_eff_line_mdot[iN][i][1] = np.nan
                        constant_eff_line_pr[iN][i][1] = np.nan
            else:
                constant_eff_line_mdot[iN][i][0] = np.nan
                constant_eff_line_pr[iN][i][0] = np.nan
                constant_eff_line_mdot[iN][i][1] = np.nan
                constant_eff_line_pr[iN][i][1] = np.nan


        # checking if all efficiencies are greater than 0 and less or equal to 1
        if any(etao < 0) or any(etao > 1):  
            print('\tError: Rotational speed  = ' + str(Narr[iN]) + ' RPM has efficiency outside range [0, 1]. Efficiencies =', etao)
            N_errors += 1
                
        Nplot = Narr[iN]
        Nplot = str(Nplot * 10 ** (- 3)) + ' krpm'
        
        # Set plot parameters------------------------------------------------------------------------------
        plt.rcParams.update(plt.rcParamsDefault)
        plt.rcParams.update({'font.size': 12})
        
        # Plotting pressure ratio 
        fig1 = plt.figure('Pressure ratio 1')
        plt.plot(mdoto, Pro, label = Nplot) # marker = markers[iN],
        plt.ylabel(r'$\frac{P_{03}}{P_{00}}$', rotation = 45, fontsize = 12)
        plt.xlabel(r'$\dot{m}$' + ' ' + '[kg/s]', fontsize = 12)      
        plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        plt.tight_layout(rect = [0, 0, 0.98, 1])
        plt.plot([0, 50], [1, 1], 'r--')
        plt.title(r'Pressure ratio 1')
        plt.grid(True)
    
        # Compressor map 
        fig99 = plt.figure('Pressure ratio 2')
        plt.plot(correctedMass, Pro, label = Nplot) # marker = markers[iN]
        if iN == 5:
            arg = np.argmax(Pro <= 1.24)
            T2obs_design_point = (T2oabs[arg] + T2oabs[arg - 1]) / 2
            plt.plot(InletConditions.mdot * Compressor.Pr / np.sqrt(T2obs_design_point / InletConditions.T00), Compressor.Pr, 'ro', label = 'Design point')
        plt.ylabel(r'$\frac{P_{03}}{P_{00}}$', rotation = 45, fontsize = 12)
        plt.xlabel(r'$\dot{m} \frac{P_{03}}{P_{00}} \sqrt{\frac{T_{00}}{T_{03}} }  $' + ' ' + '[kg/s]', fontsize = 12)
        plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        plt.title('Pressure ratio 2')
        plt.tight_layout(rect = [0, 0, 0.98, 1])
        # plt.xlim([0, 50])
        plt.grid(True)

        # Mach number
        fig19 = plt.figure('Impeller exit Mach number')
        plt.plot(mdoto, M2o, label = Nplot) # marker = markers[iN]
        plt.ylabel( r'${Ma}_2$', fontsize = 12)
        plt.xlabel(r'$\dot{m}$' + ' ' + '[kg/s]', fontsize = 12)
        plt.legend()
        plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        plt.title('Impeller Mach number')
        plt.tight_layout(rect = [0, 0, 0.98, 1])
        plt.grid(True)
        
        # Plotting efficiency curve 
        fig2 = plt.figure('Efficiency')
        plt.plot(mdoto, etao, label = Nplot) # marker = markers[iN],
        plt.grid()
        plt.ylabel(r'$\eta$', rotation = 45, fontsize = 12)
        plt.xlabel(r'$\dot{m} $'+ ' ' + '[kg/s]', fontsize = 12)
        plt.legend()
        plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        plt.tight_layout(rect = [0, 0, 0.98, 1])
        plt.title(r'Efficiency ')
        plt.grid(True)
        
        fig4 = plt.figure('Exit Tip velocity')
        plt.plot(Narr[iN], U2o, 'ko', label = Nplot) # marker = markers[iN],
        plt.grid()
        plt.xlabel(r'$N $  [rpm]', fontsize = 12)
        plt.ylabel(r'$U_2 $ [m/s]', fontsize = 12)
        plt.title('Impeller exit tip speed')
        plt.grid(True)

        fig100 = plt.figure('Performance map')
        plt.plot(mdoto, Pro, label = 'RPM = ' + str(Narr[iN]))

    fig1 = plt.figure('Pressure ratio 1', label = 'Design pt.')     # MSG: Add this to legend
    plt.plot(InletConditions.mdot, Compressor.Pr, 'ro')
    print('\nOff-design_performance calculated and plotted successfully.')
    print('\nNumber of errors in efficiency: ' + str(N_errors) + ' out of ' + str(len(Narr)) + ' rotational speeds')

    fig100 = plt.figure('Performance map')
    for i in range(len(constant_eff_lines)):
        line_1_mdot = [constant_eff_line_mdot[j][i][0] for j in range(len(constant_eff_line_mdot))]
        line_2_mdot = [constant_eff_line_mdot[j][i][1] for j in range(len(constant_eff_line_mdot))]
        line_1_pr = [constant_eff_line_pr[j][i][0] for j in range(len(constant_eff_line_pr))]
        line_2_pr = [constant_eff_line_pr[j][i][1] for j in range(len(constant_eff_line_pr))]
        plt.plot(line_1_mdot, line_1_pr, linestyle = '--', color = colors[i], label = 'Efficiency = ' + str(constant_eff_lines[i]))
        plt.plot(line_2_mdot, line_2_pr, linestyle = '--', color = colors[i])
    plt.xlabel(r'$\dot{m}$' + ' ' + '[kg/s]', fontsize = 12)      
    plt.ylabel(r'$\frac{P_{03}}{P_{00}}$', rotation = 45, fontsize = 12)
    plt.title('Performance map')
    plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
    plt.tight_layout(rect = [0, 0, 0.98, 1])
    plt.grid(True)
    #plt.show()