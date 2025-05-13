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
    
    Cm1test = 258.71 * 0.3048
    T1test = InletConditions.T00i - ((Cm1test) ** 2) / (2 * Fluid.Cp)
    print('T1test', (T1test[0] - 273.15) * 9 / 5 + 32)
    P1test = InletConditions.P00 * (T1test / InletConditions.T00) ** (Fluid.k / (Fluid.k - 1))
    print('P1test', P1test[0] * 0.000145037738)
    rho1test = P1test / (Fluid.R * T1test)  
    A1test = InletConditions.mdot / (rho1test * Cm1test * (1 - InletConditions.B1))
    print('A1test', A1test[0] / 0.0254 ** 2)
    rt1test = (A1test / (np.pi * (1 - Compressor.rhDivr1 ** 2))) ** 0.5
    print('rt1test', rt1test[0] / 0.0254)
    U1ttest = 2 * np.pi * rt1test * N / 60                                                # Inducer blade tip speed [m/s]
    print('U1ttest', U1ttest[0] / 0.3048)
    W1ttest = (Cm1test ** 2 + (U1ttest - Ctheta1i) ** 2) ** 0.5              # Inducer relative velocity [m/s]
    print('W1test', W1ttest[0] / 0.3048)

    M1i = C1i / ((Fluid.k * Fluid.R * T1i) ** 0.5)                                  # Inducer Mach number [-]                                 by definition
    P1i = InletConditions.P00 * (T1i / InletConditions.T00) ** (Fluid.k / (Fluid.k - 1))    # Inducer pressure [Pa]                           Isentropic relation                        
    rho1i = P1i / (Fluid.R * T1i)                                                   # Inducer density  [kg/m^3]                               Assuming ideal gas                          
    A1i = InletConditions.mdot / (rho1i * InletConditions.Cm1i * (1 - InletConditions.B1))  # Inducer flow area [m^2]                         Continuity                                                        
    if Compressor.optimize_inducer_geometry:        # If rh/r1 is given
        rt1i = (A1i / (np.pi * (1 - Compressor.rhDivr1 ** 2))) ** 0.5                   # Inducer tip radius [m]                                  Geometry
        rhi = rt1i * Compressor.rhDivr1
    else:   # If rh is given
        rt1i = (A1i / np.pi + Compressor.rh ** 2) ** 0.5
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
    
    print('Cm1', Compressor.Cm1 / 0.3048, 'W1t', Compressor.W1t / 0.3048, 'U1t', Compressor.U1t / 0.3048, 'P1', Compressor.P1 * 0.000145037738, 'T1', Compressor.T1, 'beta1', Compressor.beta1, 'A1', Compressor.A1 / 0.0254 ** 2, 'r1t', Compressor.r1 / 0.0254, Compressor.M1)
    #print(1/0)

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
            print('\nBlade number:', IterationMatrix.ZBarr[iz], 'Blade angle:', np.rad2deg(IterationMatrix.beta2BArr[ib][0]))
            sigma = IterationMatrix.sigmaWiesnerMat[ib, iz]     # slip factor for given blade number and blade angle
            etaStage = Compressor.etaStage0                     # resetting the guessed efficiency for each new combination of blade number and blade angle
            Pr_stage = 0
            while Pr_stage / Compressor.Pr < 0.999 or Pr_stage / Compressor.Pr > 1.001:
                dh0s = ((Fluid.k * Fluid.R * InletConditions.T00) / (Fluid.k - 1)) * (Compressor.Pr ** ((Fluid.k - 1) / Fluid.k) - 1)       
                Wx = dh0s / etaStage                                                # Specific work [J/kg/K]
                
                T02m = InletConditions.T00 + Wx * (Fluid.k - 1) / (Fluid.k * Fluid.R)        # Stagnation exit temperature [K]     , from dh=cp*Dt   
                mu = sigma * Compressor.lambda2 / (Compressor.lambda2 - np.tan((IterationMatrix.beta2BArr[ib][0])))
                U2 = ((Compressor.U1t * Compressor.Ctheta1 + Wx) / mu) ** (1 / 2)   # Impeller tip velocity [m/s]    , from work input coefficient MSG: Should be divided by mdot?
                D2 = 60 * U2 / (np.pi * Compressor.N0)
                r2 = D2 / 2
                Ctheta2m = mu * U2
                Cm2m = Ctheta2m / Compressor.lambda2
                C2 = (Ctheta2m ** 2 + Cm2m ** 2) ** 0.5 
                T2m = T02m - (Fluid.k - 1) / (2 * Fluid.k * Fluid.R) * C2 ** 2      # Exit temperature [K]                            from stagnation temperature            
                M2 = U2 / np.sqrt(Fluid.k * Fluid.R * T2m)                          # Impeller exit blade mach number                 from definition
                P02m = InletConditions.P00 * ((Wx * Compressor.eta_rotor * (Fluid.k - 1)  / (Fluid.k * Fluid.R * InletConditions.T00)) + 1) ** (Fluid.k / (Fluid.k - 1))     # Exit Stagnation Pressure, isentropic [Pa]                   
                P2m = P02m / ((T02m / T2m) ** (Fluid.k / (Fluid.k - 1)))            # Exit pressure [Pa]                              from stagnation pressure, isentropic relation
                rho2m = P2m / (T2m * Fluid.R)                                       # Exit density [kg/m^3]                           from ideal gas law
                A2m = InletConditions.mdot / (rho2m * Cm2m)                         # Exit area [m^2]                                 from continuity
                b2 = A2m / (np.pi * D2)                                             # Impeller exit cylinder height [m]               from geometry
                
                CpDi = 1 - 1 / Compressor.AR ** 2   
                CpD = Compressor.etad * CpDi
                
                p5 = P2m + CpD * (P02m - P2m)                                       # From Bernoulli's eq. for a vaneless diffuser

                etaStage = ((p5 / InletConditions.P00) ** ((Fluid.k - 1) / Fluid.k) - 1) / (T02m / InletConditions.T00 - 1)
                Pr_stage = p5 / InletConditions.P00
                print('etaStage', etaStage)
                print('Pr_stage', Pr_stage)
                print('T02m', T02m)
                
        


            IterationMatrix.etaMat[ib, iz] = etaStage
            IterationMatrix.MachExitMat[ib, iz] = M2
            IterationMatrix.b2Mat[ib, iz] = b2
            IterationMatrix.PrestMat[ib, iz] = p5 / InletConditions.P00
            IterationMatrix.WxMat[ib, iz] = Wx
            IterationMatrix.c2Mat[ib, iz] = C2
            IterationMatrix.c2mMat[ib, iz] = Cm2m
            IterationMatrix.Ctheta2Mat[ib, iz] = Ctheta2m
            IterationMatrix.sigmaMat[ib, iz] = sigma
            IterationMatrix.r2Mat[ib, iz] = r2  
            IterationMatrix.U2Mat[ib, iz] = U2
            print('first', (Fluid.k - 1) / (2 * Fluid.k * Fluid.R), 'second', C2 ** 2, 'D2', D2, 'U1t', Compressor.U1t, 'rh1', Compressor.rh1, 'r1', Compressor.r1, 'rh/r1', Compressor.rh1 / Compressor.r1)

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
    else:
        Compressor.r2 = IterationMatrix.r2Mat[np.where(IterationMatrix.beta2BArr == Compressor.bladeAngle)[0][0]][np.where(IterationMatrix.ZBarr == Compressor.bladeNumber)[0][0]]
        #print('T02m', (IterationMatrix.U2Mat[np.where(IterationMatrix.beta2BArr == Compressor.bladeAngle)[0][0]][np.where(IterationMatrix.ZBarr == Compressor.bladeNumber)[0][0]] - 273.15) * 9 / 5 + 32)   
        #print('U2', IterationMatrix.U2Mat[np.where(IterationMatrix.beta2BArr == Compressor.bladeAngle)[0][0]][np.where(IterationMatrix.ZBarr == Compressor.bladeNumber)[0][0]] / 0.3048)
        
        print('T2', IterationMatrix.etaMat[np.where(IterationMatrix.beta2BArr == Compressor.bladeAngle)[0][0]][np.where(IterationMatrix.ZBarr == Compressor.bladeNumber)[0][0]])
        Compressor.D2 = 2 * Compressor.r2 
        Compressor.U2 = IterationMatrix.U2Mat[np.where(IterationMatrix.beta2BArr == Compressor.bladeAngle)[0][0]][np.where(IterationMatrix.ZBarr == Compressor.bladeNumber)[0][0]]
        Compressor.b2 = IterationMatrix.b2Mat[np.where(IterationMatrix.beta2BArr == Compressor.bladeAngle)[0][0]][np.where(IterationMatrix.ZBarr == Compressor.bladeNumber)[0][0]]
        #print(1/0)
        print('\nGeometry successfully calculated')                


def print_and_plot_geometry(Compressor, InletConditions, IterationMatrix):
    # ------------- Preparing data for plotting ------------ 
    """ Plotting functions plotCompressor.py, plotFlow.py, plotSystem.py, plotVelocities.py and plotText.py are 
            called from logic.py but prepared here since data is produced here."""

    #print('Critical RPM: ' + str(round(Compressor.Ncrit, 2)) )
    print('Applied RPM: ' + str(round(Compressor.Ndes, 2)) )
    print('rh: ' + str(round(Compressor.rh1, 10)) + 'm' )
    print('r1: ' + str(round(Compressor.r1, 10)) + 'm')
    print('r2: ' + str(round(Compressor.r2 , 10)) + 'm' )
    print('b2: ' + str(round(Compressor.b2 , 10)) + 'm' )

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
            r"$U_{2t,crit}= $" + str(round(Compressor.U2Crit, 3)) +"m/s"+ "\n" 
    text2 = "\n" \
            r"$r_{t1}$ = " + str(round(Compressor.r1, 3)) + r"$m$" +" \n" \
            r"$r_{h1}$ = " + str(round(Compressor.rh1, 3)) + r"$m$" +" \n" \
            r"$r_{t2}$ = " + str(round(Compressor.r2, 3)) + r"$m$" +" \n" \
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

if __name__ == '__main__':
    pass