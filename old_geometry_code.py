"""
The following Python code takes a set of variables as inputs and calculates the geometry of an inducer, 
impeller and diffuser for a centrifugal compressor.

Reference: Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244
Author(s): Petter Resell (summer intern, 2024), Martin Spillum Grønli (SINTEF Energy Research, 2024)

"""

### Import-------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

from pressure_test import pressureOverUnderEstimate
import geometry_system_functions


def minimize_relative_velocity(N, Fluid, InletConditions, Compressor):
    """ Function to find the inducer properties which minimizes the relative velocity for a given rotational speed N. """
    
    # Make arrays for all possible values of Cm1
    Ctheta1i = InletConditions.Cm1i * np.tan(np.radians(InletConditions.alpha1))    # Inducer absolute tangential velocity [degrees]          Trigonometry
    C1i = (Ctheta1i ** 2 + InletConditions.Cm1i ** 2) ** 0.5                        # Inducer Absolute velocity C1 [m/s]                      Pythagoras theorem from velocity triangle 
    T1i = InletConditions.T00i - (C1i ** 2) / (2 * Fluid.Cp)                        # Inducer temperature [K]                                 Stagnation temperature relation             
    M1i = C1i / ((Fluid.k * Fluid.R * T1i) ** 0.5)                                  # Inducer Mach number [-]                                 by definition
    P1i = InletConditions.P00 * (T1i / InletConditions.T00) ** (Fluid.k / (Fluid.k - 1))    # Inducer pressure [Pa]                           Isentropic relation                        
    rho1i = P1i / (Fluid.R * T1i)                                                   # Inducer density  [kg/m^3]                               Assuming ideal gas                          
    A1i = InletConditions.mdot / (rho1i * InletConditions.Cm1i * (1 - InletConditions.B1))  # Inducer flow area [m^2]                         Continuity                                                        
    rt1i = (A1i / (np.pi * (1 - Compressor.rhDivr1 ** 2))) ** 0.5                   # Inducer tip radius [m]                                  Geometry
    U1ti = 2 * np.pi * rt1i * N / 60                                                # Inducer blade tip speed [m/s]
    W1ti = (InletConditions.Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5              # Inducer relative velocity [m/s]
    InletConditions.W1ti = W1ti                                                     # Store relative velocities for plotting

    # Get the index of the minimum relative velocity
    w1ti_index_min = np.argmin(W1ti)                                                # Index of smallest relative velocity [m/s]
    
    # Set inducer properties to the values that minimize the relative velocity
    Compressor.r1 = rt1i[w1ti_index_min]                                            # Inducer tip radius [m]
    Compressor.Ctheta1 = Ctheta1i[w1ti_index_min]                                   # Inducer angular velocity [m/s]
    Compressor.C1 = C1i[w1ti_index_min]                                             # Inducer velocity [m/s]
    Compressor.T1 = T1i[w1ti_index_min]                                             # Inducer temperature [K]
    Compressor.M1 = M1i[w1ti_index_min]                                             # Inducer Mach number [-]
    Compressor.P1 = P1i[w1ti_index_min]                                             # Inducer pressure [Pa]
    Compressor.rho1 = rho1i[w1ti_index_min]                                         # Inducer density [kg/m3]
    Compressor.A1 = A1i[w1ti_index_min]                                             # Inducer flow area [m2]
    Compressor.U1t = U1ti[w1ti_index_min]                                           # Inducer tip speed [m/s]
    Compressor.Cm1 = InletConditions.Cm1i[w1ti_index_min]                           # Inducer meridional velocity [m/s]
    Compressor.beta1 = np.degrees(np.arctan((Compressor.U1t - Compressor.Ctheta1) / Compressor.Cm1))    # Inducer relative velocity angle [deg]
    Compressor.W1t = W1ti[w1ti_index_min]                                           # Inducer relative velocity [m/s]
    Compressor.omega = Compressor.U1t / Compressor.r1                               # Angular velocity [rad/s]   


def inducer_and_impeller_calculations(Fluid, InletConditions, Compressor):
    """ Impeller calculations """
    
    """ The isentropic enthalpy demand is used to find the approximate required work. This is done through iterating the isentropic efficiency.
        Furthermore the  Work is used to find the outlet stagnation temperature which in turn is used to find the outlet stagnation pressure 
        from isentropic relations. Furthermore, geometrical relations are applied. """

    """ If the pressure estimate Prest initially over-estimates the pressure ratio by more than 5%, then the the error will only grow for increasing efficiency
            since the pressure ratio increase for increasing efficiency. Hence its necessary to check if the pressure estimate is an 
                under- or over-estimate before adjusting the efficiency eta. From then the efficiency can be either increased or decreased. """

    """ Iterating to find a impeller tip velocity that satisfy material constraint. Taking the rotational speed set in settings.py and
            checking if the corresponding rotational velocity U2 is bigger than the critical rotational velocity. If it is then 
                inducerFlowWRTminimumRelativeVelocity(N=Ndes) finds all inducer flow properties for the new N and U2. For the first iteration U2=U_crit 
                so it runs at least one iteration anyways. If the new U2 is found to be to large then it is lowered by increments of 1000.
                    This could be replaced by seting the rotational velocity to for example 70% or 80% of the critical one.  """

    Compressor.Ndes = Compressor.N0        # Initial guess for N    
    minimize_relative_velocity(N = Compressor.Ndes, Fluid = Fluid, InletConditions = InletConditions, Compressor = Compressor)
    Compressor.rh1 = Compressor.rhDivr1 * Compressor.r1                         # Hub radius [m]
    
    """
    while (U2) >= Compressor.U2Crit and Ndes > 0:       # MSG: Find a way to find the lowest rotational speed, instead of the highest?    

        if (Compressor.U2Crit) > Compressor.bladeVelUpperLimit:      # Can skip this break to avoid crashing code
            break
        if (U2) >= Compressor.U2Crit:
            Ndes -= 1000
    """
    

def checkUpdateSlipFactorSigma(epsilonLimit, slipFactor, Compressor):
    """ Function to check if condition for the use of wiesner slip factor is upheld. If its not upheld then slip factor 
            is updated by the correction equation. If the condition is upheld it just returns the same slip factor. 
            The Function is implemented when iterating through all blade numbers and blade angles below.  """
    if Compressor.r1Divr2 >= epsilonLimit:
        return slipFactor * (1 - ((Compressor.r1Divr2 - epsilonLimit) / (1 - epsilonLimit)) ** 3)
    else:
        return slipFactor


def iterate_blade_number_and_blade_angle(Compressor, InletConditions, Fluid, IterationMatrix):
    """ Iterating main functionality of design procedure. Looping through each blade number and blade angle to 
            capture all combinations.  """    

    for iz in range(len(IterationMatrix.ZBarr)):
        for ib in range(len(IterationMatrix.beta2BArr)):
            #print('Blade number:', IterationMatrix.ZBarr[iz], 'Blade angle:', np.rad2deg(IterationMatrix.beta2BArr[ib][0]))
            sigma = IterationMatrix.sigmaWiesnerMat[ib, iz]     # slip factor for given blade number and blade angle
            etaStage = Compressor.etaStage0                     # resetting the guessed efficiency for each new combination of blade number and blade angle
            trueFalseCheck = False                              # iteration control for while loop
            
            """ Handling slip factor. Updating if the conditions for wiesner relation is not upheld. """
            epsilonLimit = IterationMatrix.rhsExpLimitMat[ib, iz]                                                           # right hand side of wiesner slip factor condition
            sigma = checkUpdateSlipFactorSigma(epsilonLimit = epsilonLimit, slipFactor = sigma, Compressor = Compressor)    # updating slip factor wrt wiesner condition. MSG: Comment out?
            
            """ impellerOutletVelocities() finds the impeller outlet velocities and work output coeff. """
            #workInputCoeff, Ctheta2m, Cm2m, C2 = geometry_system_functions.impellerOutletVelocities(slipFactor = sigma, beta2B = (IterationMatrix.beta2BArr[ib])[0], U2 = Compressor.U2, lambda2 = Compressor.lambda2)          # Finding impeller outlet velocities     # MSG: Function need to be updated with classes

            """ The following five lines are made to find porperties of the slip for a given blade number and blade angle.
                    None of the resulting variables are applied in further calculations, but could easily be. They are included
                    for further development and demonstration of slip.  """
            
            #Ctheta2ideal = Compressor.U2 - Ctheta2m * np.tan((np.abs((IterationMatrix.beta2BArr[ib])[0])))   
            #CTheta2Real = sigma * Ctheta2ideal         
            #Cslip1 = Ctheta2ideal - CTheta2Real
            #beta2flow = np.rad2deg(np.arctan((Cslip1 + Cm2m * np.tan(np.abs(IterationMatrix.beta2BArr[ib]))) / Cm2m))
            #dh0SlipCorrected =  Compressor.U2 * (CTheta2Real) - Compressor.U1t * Compressor.Ctheta1   # Alternative for finding fluid enthalpy change, for comparison

            """ The while loop iterates for a given combination of blade number and blade angle. It iterates the efficiency demand
                    while the estimated pressure ratio is not within the given tolerance of the desired pressure ratio. It is updated in 
                    the function called pressureOverUnderEstimate. The guess efficiency is increased or decreased varying if the 
                    pressure ratio is under- or overestimated. """
            Pr_goal = Compressor.Pr
            Pr_isentropic = Compressor.Pr
            Pr = Compressor.Pr
            stopIteration = False
            eta_max = 0.99999
            while stopIteration == False:# and Compressor.U2 < Compressor.bladeVelUpperLimit:
                # -------Fluid------------ Updating work with the new isentropic efficiency -------------------
                print(iz)
                print('Pr_goal:', Pr_goal)
                print('etaStage:', etaStage)
                Compressor.dh0s = ((Fluid.k * Fluid.R * InletConditions.T00) / (Fluid.k - 1)) * ((Pr_isentropic) ** ((Fluid.k - 1) / Fluid.k) - 1)       
                Wx = Compressor.dh0s / etaStage                         # Specific work [J/kg/K]
                #workError = np.abs(Wx - dh0SlipCorrected) / Wx          # Comparing the two methods of finding enthalpy change
                # Wx = dh0SlipCorrected                                   # Work set to be given by slip corrected euler equation                 

                """ Impeller outlet calculation. Stagnation properties denoted by zero. Isentropic relations applied on next four rows. """
                T02m = InletConditions.T00 + Wx * (Fluid.k - 1) / (Fluid.k * Fluid.R)        # Stagnation exit temperature [K]     , from dh=cp*Dt   
                mu = sigma * Compressor.lambda2 / (Compressor.lambda2 - np.tan(np.abs(IterationMatrix.beta2BArr[ib])))
                Compressor.U2 = ((Compressor.U1t * Compressor.Ctheta1 + Wx) / (mu)) ** (1 / 2)    # Impeller tip velocity [m/s]    , from work input coefficient MSG: Should be divided by mdot?
                Compressor.D2 = 60 * Compressor.U2 / (np.pi * Compressor.N0)
                print('N, U2, D2', Compressor.N0, Compressor.U2, Compressor.D2)
                Compressor.r2 = Compressor.D2 / 2
                print('rh1', Compressor.rh1, 'r1', Compressor.r1, 'r2', Compressor.r2)
                print('rh1/r1', Compressor.rh1 / Compressor.r1, 'r1/r2', Compressor.r1 / Compressor.r2)
                Ctheta2m = mu * Compressor.U2
                Cm2m = Ctheta2m / Compressor.lambda2
                C2 = (Ctheta2m ** 2 + Cm2m ** 2) ** 0.5 
                T2m = T02m - (Fluid.k - 1) / (2 * Fluid.k * Fluid.R) * C2 ** 2             # Exit temperature [K]                            from stagnation temperature            
                P02m = InletConditions.P00 * ((Compressor.dh0s * (Fluid.k - 1)  / (Fluid.k * Fluid.R * InletConditions.T00)) + 1) ** (Fluid.k / (Fluid.k - 1))     # Exit Stagnation Pressure, isentropic [Pa]                   
                P2m = P02m / ((T02m / T2m) ** (Fluid.k / (Fluid.k - 1)))                     # Exit pressure [Pa]                              from stagnation pressure, isentropic relation
                rho2m = P2m / (T2m * Fluid.R)                                                # Exit density [kg/m^3]                           from ideal gas law
                A2m = InletConditions.mdot / (rho2m * Cm2m)                                   # Exit area [m^2]                                 from continuity
                b2 = A2m / (np.pi * Compressor.D2)                                            # Impeller exit cylinder height [m]               from geometry
                
                M2 = Compressor.U2 / np.sqrt(Fluid.k * Fluid.R * T02m)                       # Impeller exit blade mach number                 from definition

                # ------------------- Finding diffuser properties -------------------
                P3, P03, C3 = geometry_system_functions.diffuserFlow(P2 = P2m, P02 = P02m, rho2 = rho2m, C2 = C2, CpD = Compressor.CpD, AR = Compressor.AR)   # Finding diffuser properties
                Pr = P3 / InletConditions.P00
                etaCalculated = ((P03 / InletConditions.P00) ** ((Fluid.k - 1) / Fluid.k) - 1) / ((T02m / InletConditions.T00) - 1)        # Iterative stage efficiency [-], isothermal diffuser assumed?
                Prest = ((etaCalculated * Compressor.U2 ** 2 * mu) / (Fluid.Cp * Compressor.T1) + 1) ** (Fluid.k / (Fluid.k - 1))          # Estimate of the pressure ratio, equation is validated
                print('Pr:', Pr, 'Prest', Prest, 'etaStage', etaStage, 'etaCalculated', etaCalculated)

                #etaStage = etaCalculated
                #if etaCalculated > Compressor.etaUpperLimit or etaCalculated < Compressor.etaLowerLimit:  # MSG: Do this after the pressure test? This will probably stop during the first run since etaCalculated will too high
                #    stopIteration = True
                #    print('etaCalculated:', etaCalculated)
                #    break
                #else:
                if np.abs((P3 / InletConditions.P00) - Pr_goal) < Compressor.iterTol:
                    print('yes')
                    stopIteration = True
                    IterationMatrix.etaMat[ib, iz] = etaStage#etaCalculated
                    #IterationMatrix.pressErrorMat[ib, iz] = abs(PressureTestOuterLoop)
                    IterationMatrix.MachExitMat[ib, iz] = M2
                    IterationMatrix.b2Mat[ib, iz] = b2
                    IterationMatrix.PrestMat[ib, iz] = Pr#Prest
                    IterationMatrix.WxMat[ib, iz] = Wx
                    #IterationMatrix.dh0SlipCorrMAt[ib, iz] = dh0SlipCorrected
                    IterationMatrix.c2Mat[ib, iz] = C2
                    IterationMatrix.c2mMat[ib, iz] = Cm2m
                    IterationMatrix.Ctheta2Mat[ib, iz] = Ctheta2m
                    #IterationMatrix.VslipMat[ib, iz] = Cslip1
                    IterationMatrix.sigmaMat[ib, iz] = sigma
                    #IterationMatrix.beta2flowMat[ib, iz] = beta2flow
                    #continue
                elif ((P3 / InletConditions.P00) - Pr_goal) < 0:
                    etaStage += 0.005
                    if etaStage > eta_max:
                        Pr_goal -= 0.1
                        Pr_isentropic += 0.1
                        etaStage = Compressor.etaStage0
                        if Pr_goal <= 1:
                            stopIteration = True
                elif ((P3 / InletConditions.P00) - Pr_goal) > 0:
                    etaStage -= 0.005
                
                
                    
                    

                """ Updating efficiency for next iteration"""
                #etaStage = pressureOverUnderEstimate(PressureTestOuterLoop, etaStage)   # MSG: Function needs to be updated with classes
        

    # Check if any valid geometries were found
    countTrue = np.count_nonzero(~np.isnan(IterationMatrix.etaMat))  
    if countTrue == 0:
        # If no valid geometries were found for the pressure ratio, lower the pressure ratio and try again
        IterationMatrix.etaMat[:] = np.nan          # Reset all values to nan
        IterationMatrix.pressErrorMat[:] = np.nan
        IterationMatrix.MachExitMat[:] = np.nan
        IterationMatrix.b2Mat[:] = np.nan
        IterationMatrix.PrestMat[:] = np.nan
        IterationMatrix.WxMat[:] = np.nan
        IterationMatrix.dh0SlipCorrMAt[:] = np.nan
        IterationMatrix.c2Mat[:] = np.nan
        IterationMatrix.c2mMat[:] = np.nan
        IterationMatrix.Ctheta2Mat[:] = np.nan
        IterationMatrix.VslipMat[:] = np.nan
        IterationMatrix.sigmaMat[:] = np.nan
        IterationMatrix.beta2flowMat[:] = np.nan

        if Compressor.Pr <= 1.0:
            print(f"Combination of rh/r1, r1/r2 and N gave zero valid cases for pressure ratios greater than 1.0")
            #raise Exception(f"Combination of rh/r1, r1/r2 and N gave zero valid cases for pressure ratios greater than 1.0")
    else:
        print('\nGeometry successfully calculated')                


def print_and_plot_geometry(Compressor, InletConditions, IterationMatrix):
    # ------------- Preparing data for plotting ------------ 
    """ Plotting functions plotCompressor.py, plotFlow.py, plotSystem.py, plotVelocities.py and plotText.py are 
            called from logic.py but prepared here since data is produced here."""

    #print('Critical RPM: ' + str(round(Compressor.Ncrit, 2)) )
    print('Applied RPM: ' + str(round(Compressor.Ndes, 2)) )
    print('Rotational velocity: ' + str(round(Compressor.U2[0], 2)))
    print('rh: ' + str(round(Compressor.rh1, 10)) + 'm' )
    print('r1: ' + str(round(Compressor.r1, 10)) + 'm')
    print('r2: ' + str(round(Compressor.r2[0], 10)) + 'm' )

    print('\nPlotting geometry...')

    """ Making text to go with plots """
    maxEta = np.nanmax(IterationMatrix.etaMat)
    minEta = np.nanmin(IterationMatrix.etaMat)
    countTrue = np.count_nonzero(~np.isnan(IterationMatrix.etaMat))   # MSG: Understand this condition
    totalCases = len(IterationMatrix.ZBarr) * len(IterationMatrix.beta2BArr)
    print('Valid cases: ' + str(countTrue) + '/'  +str(totalCases))
    text1 = "\n" \
            r"Valid cases: " + str(countTrue) + '/'  +str(totalCases) + "\n" \
            r"$\eta _{max}$ = " + str(round(maxEta, 3)) + "   " \
            r"$\eta _{min}$ = " + str(round(minEta, 3)) + "\n" \
            r"$W_{1t}= $" + str(round(Compressor.W1t, 3)) + "m/s   " \
            r"$C_{m1}= $" + str(round(Compressor.Cm1, 3)) + "m/s"+"\n" \
            r"$U_{1t}= $" + str(round(Compressor.U1t, 3)) + "m/s"+"   " \
            r"$M_{1f}= $" + str(round(Compressor.M1, 3)) + "\n" \
            r"$U_{2t}$ = "  + str(round(Compressor.U2[0], 3)) + "m/s"+"   " \
            r"$U_{2t,crit}= $" + str(round(Compressor.U2Crit, 3)) +"m/s"+ "\n" 
    text2 = "\n" \
            r"$r_{t1}$ = " + str(round(Compressor.r1, 3)) + r"$m$" +" \n" \
            r"$r_{h1}$ = " + str(round(Compressor.rh1, 3)) + r"$m$" +" \n" \
            r"$r_{t2}$ = " + str(round(Compressor.r2[0], 3)) + r"$m$" +" \n" \
            #r"$N_{crit}$ = " + str(round(Compressor.Ncrit, 3)) + r"$rpm$" +" \n"  
    text3 = "\n" \
            r"Desired $PR^*$: " + str(Compressor.Pr) +  "\n" \
            r"$\frac{r_h}{r_1}{*}$ = " + str(Compressor.rhDivr1) + r" $\frac{m}{m} $" +"    " \
            r"$\frac{r_1}{r_2}{*}$ = " + str(Compressor.r1Divr2) + r" $\frac{m}{m} $" +"\n" \
            r"$N^{*}$ = " +str(round(Compressor.Ndes, 1)) + r"$rpm$" +" \n" \
            r"Proposed efficiency $\eta *=$  " + str(Compressor.etaStage0) + "\n" \
            r"Exit swirl number $\lambda_2 *$: " + str(Compressor.lambda2) +  "\n" \
            r"Tolerance*: " + str(round(Compressor.iterTol, 3)) + "   "
    text = [text1, text2, text3]     # Used in plotText.py

    """ Plot of relative velocities as a function of inlet meridional velocity """
    fig, ax11 = plt.subplots()
    plt.plot(InletConditions.Cm1i, InletConditions.W1ti)
    ax11.tick_params(axis = 'both', which = 'major', labelsize = 10)
    plt.xlabel(r'$C_{\mathrm{m1}}$ [m/s]', fontsize = 12)
    plt.ylabel(r'$W_{\mathrm{1t}}$ [m/s]', fontsize = 12)
    plt.title(r'Minimization of relative velocity $W_{\mathrm{1t}}$')
    plt.grid()

    # Import plotting scripts
    from plot_scripts.plot_compressor import plotCompressorParam
    from plot_scripts.plot_flow import plotFlowConditions
    from plot_scripts.plot_velocities import plotVelocities
    from plot_scripts.plot_system import plotSystemVariables
    from plot_scripts.plot_text import plotText
    
    # Set plot parameters------------------------------------------------------------------------------
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'font.size': 15})
    #plt.close('all')
    
    # ------------------- Plotting from geometry.py -------------------
    plotCompressorParam(Compressor = Compressor, IterationMatrix = IterationMatrix)
    plotFlowConditions(Compressor = Compressor, IterationMatrix = IterationMatrix)
    plotVelocities(Compressor = Compressor, IterationMatrix = IterationMatrix)
    plotSystemVariables(Compressor = Compressor, IterationMatrix = IterationMatrix)
    plotText(text)
    #plt.show()