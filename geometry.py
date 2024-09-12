"""
The following Python code takes a set of variables as inputs and calculates the geometry of an inducer, 
impeller and diffuser for a centrifugal compressor.

Reference: Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244
Author: Petter Resell (summer intern, 2024)

"""

### Import-------------------------------------------------------------------------------------------
import math
import numpy as np
import matplotlib.pyplot as plt


from pressure_test import pressureOverUnderEstimate
import geometry_system_functions

plt.rcParams.update(plt.rcParamsDefault)            # For plotting
plt.rcParams.update({'font.size': 15})

print('\n' + 'Running Geometry.py')  


# Declaring variables for off-design.py. MSG: I think these can be removed or moved to off-design_performance.py
#U1t = 0     # MSG: Try deleting this. Should now be part of Compressor class
#beta1 = 0


def inducer_calculations(Fluid, InletConditions, Compressor):
    """ Inducer calculations
    For each value of Cm1 in settings.Cm1i, we calculate the corresponding values of C1, T1, M1, P1, rho1, A1, r1, U1t, W1t, Cm1 and beta1.
        The value of Cm1 will affect the values of U1t and Ctheta1 (through ex. r1), which in turn affects the value of W1t. Isentropic relations 
            linking the static to the total conditions. Velocity triangles are also applied. 
                By plotting W1t against Cm1, we then find the value of Cm1 that minimises W1t. """
    Compressor.Ctheta1i = InletConditions.Cm1i * math.tan(math.radians(InletConditions.alpha1))                         # Inducer absolute tangential velocity [degrees]          Trigonometry
    Compressor.C1i = (Compressor.Ctheta1i ** 2 + InletConditions.Cm1i ** 2) ** 0.5                                      # Inducer Absolute velocity C1 [m/s]                      Pythagoras theorem from velocity triangle 
    Compressor.T1i = InletConditions.T00i - (Compressor.C1i ** 2) / (2 * Fluid.Cp)                                      # Inducer temperature [K]                                 Stagnation temperature relation             
    Compressor.M1i = Compressor.C1i / ((Fluid.k * Fluid.R * Compressor.T1i) ** 0.5)                                     # Inducer Mach number [-]                                 by definition
    Compressor.P1i = InletConditions.P00 * (Compressor.T1i / InletConditions.T00) ** (Fluid.k / (Fluid.k - 1))          # Inducer pressure [Pa]                                   Isentropic relation                        
    Compressor.rho1i = Compressor.P1i / (Fluid.R * Compressor.T1i)                                                      # Inducer density  [kg/m^3]                               Assuming ideal gas                          
    Compressor.A1i = InletConditions.mdot / (Compressor.rho1i * InletConditions.Cm1i * (1 - InletConditions.B1))        # Inducer flow area [m^2]                                 Continuity                                                        
    Compressor.rt1i = (Compressor.A1i / (math.pi * (1 - Compressor.rhDivr1 ** 2) )) ** 0.5                              # Inducer tip radius [m]                                  Geometry


def inducerFlowWRTminimumRelativeVelocity(N, Compressor, InletConditions):
        """ This function takes a rotational speed N and findes the relative velocity. This is then minimized to increase efficiency. 
            Inducer properties and velocities are then found at the point where the relative velocity is the smallest. 
            This is at w1ti_index_min """
        Compressor.U1ti = 2 * math.pi * Compressor.rt1i * N / 60                                                # Inducer blade tip speed [m/s]
        Compressor.W1ti = (InletConditions.Cm1i ** 2 + (Compressor.U1ti - Compressor.Ctheta1i) ** 2) ** 0.5     # Inducer relative velocity [m/s]
        Compressor.w1ti_index_min = np.argmin(Compressor.W1ti)                                                  # Index of smallest relative velocity [m/s]     # MSG: I think this guarantees the optimal Cm1
        Compressor.r1 = Compressor.rt1i[Compressor.w1ti_index_min]                                              # Inducer tip radius [m]
        Compressor.Ctheta1 = Compressor.Ctheta1i[Compressor.w1ti_index_min]                                     # Inducer angular velocity [m/s]
        Compressor.C1 = Compressor.C1i[Compressor.w1ti_index_min]                                               # Inducer velocity [m/s]
        Compressor.T1 = Compressor.T1i[Compressor.w1ti_index_min]                                               # Inducer temperature [K]
        Compressor.M1 = Compressor.M1i[Compressor.w1ti_index_min]                                               # Inducer Mach number [-]
        Compressor.P1 = Compressor.P1i[Compressor.w1ti_index_min]                                               # Inducer pressure [Pa]
        Compressor.rho1 = Compressor.rho1i[Compressor.w1ti_index_min]                                           # Inducer density [kg/m3]
        Compressor.A1 = Compressor.A1i[Compressor.w1ti_index_min]                                               # Inducer flow area [m3]
        Compressor.U1t = Compressor.U1ti[Compressor.w1ti_index_min]                                             # Inducer tip speed [m/s]
        Compressor.Cm1 = InletConditions.Cm1i[Compressor.w1ti_index_min]                                        # Inducer meridonal velocity [m/s]
        Compressor.beta1 = math.degrees(math.atan((Compressor.U1t - Compressor.Ctheta1) / Compressor.Cm1))                 # Inducer relative velocity angle [deg]
        Compressor.W1t = Compressor.W1ti[Compressor.w1ti_index_min]                                             # Inducer relative velocity [m/s]
        Compressor.omega = Compressor.U1t / Compressor.r1                                                                  # Angular velocity [rad/s]          =(2*np.pi*N/60)



def impeller_calculations(Fluid, InletConditions, Compressor):
    """ Impeller calculations """
    
    """ The isentropic enthalpy demand is used to find the approximate required work. This is done through iterating the isentropic efficiency.
        Furthermore the  Work is used to find the outlet stagnation temperature which in turn is used to find the outlet stagnation pressure 
        from isentropic relations. Furthermore, geometrical relations are applied. """

    """ If the pressure estimate Prest initially over-estimates the pressure ratio by more than 5%, then the the error will only grow for increasing efficiency
            since the pressure ratio increase for increasing efficiency. Hence its necessary to check if the pressure estimate is an 
                under- or over-estimate before adjusting the efficiency eta. From then the efficiency can be either increased or decreased. """


    """ Specific isentropic compression enthalpy, constant value [kJ/kg/K]. dh0s is constant throughout the full iteration scheme. """
    Compressor.dh0s = ((Fluid.k * Fluid.R * InletConditions.T00) / (Fluid.k - 1) ) * (Compressor.Pr ** ((Fluid.k - 1) / Fluid.k) - 1)                     



    """ Calculating critical property of a rotating disk to use as constraint for RPM/rotational velocity/radiusettings. Disk will break at outermost point, therefore r2 and U2. """                                                                                                 
    U2Crit = np.sqrt(2 * Compressor.impellerTensileStrength / Compressor.impellerDensity)      # Applying tensile strength of disk. Titanium used.        
    Compressor.U2Crit = U2Crit
    U2 = U2Crit        
                                                                                    # initializing for loop 

    """ Iterating to find a impeller tip velocity that satisfy material constraint. Taking the rotational speed set in settings.py and
            checking if the corresponding rotational velocity U2 is bigger than the critical rotational velocity. If it is then 
                inducerFlowWRTminimumRelativeVelocity(N=Ndes) finds all inducer flow properties for the new N and U2. For the first iteration U2=U_crit 
                so it runs at least one iteration anyways. If the new U2 is found to be to large then it is lowered by increments of 1000.
                    This could be replaced by seting the rotational velocity to for example 70% or 80% of the critical one.  """
    
    Ndes = Compressor.N0
    while (U2) >= U2Crit and Ndes > 0:      # MSG: Does this give the optimal Cm1? I rememeber that this was plotted.

        inducerFlowWRTminimumRelativeVelocity(N = Ndes, Compressor = Compressor, InletConditions = InletConditions) 
        U2 = (Compressor.U1t/Compressor.r1Divr2)

        if (U2Crit) > Compressor.bladeVelUpperLimit:      # Can skip this break to avoid crashing code
            break
        if (U2) >= U2Crit:
            Ndes -= 1000
    Compressor.U2 = U2
    Compressor.Ndes = Ndes      # MSG: Change this to N0 instead of Ndes?

    """ Can find design parameters after finding U2"""
    Compressor.r2 = Compressor.U2 / Compressor.omega                          # Impeller tip radius
    Compressor.D2 = 2 * Compressor.r2                              # Impeller tip diameter
    Compressor.rh1 = Compressor.rhDivr1 * Compressor.r1            # Hub radius
    Compressor.Ncrit = 60 * U2Crit / (2 * np.pi * Compressor.r2)   # Critical rpm for the newly found impeller tip radius. 


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
        # debug=1         # breakpoint for debugging
        
        for ib in range(len(IterationMatrix.beta2BArr)):

            sigma = IterationMatrix.sigmaWiesnerMat[ib, iz]    # slip factor for given blade number and blade angle
            etaStage = Compressor.etaStage0                     # resetting the guessed efficiency for each new combination of blade number and blade angle
            trueFalseCheck = False                              # iteration control for while loop
            

            """ Handling slip factor. Updating if the conditions for wiesner relation is not upheld. """
            epsilonLimit = IterationMatrix.rhsExpLimitMat[ib, iz]                                                          # right hand side of wiesner slip factor condition
            sigma = checkUpdateSlipFactorSigma(epsilonLimit = epsilonLimit, slipFactor = sigma, Compressor = Compressor)            # updating slip factor wrt wiesner condition
            
            """ impellerOutletVelocities() finds the impeller outlet velocities and work output coeff. """
            workInputCoeff, Ctheta2m, Cm2m, C2 = geometry_system_functions.impellerOutletVelocities(slipFactor = sigma, beta2B = (IterationMatrix.beta2BArr[ib])[0], U2 = Compressor.U2, lambda2 = Compressor.lambda2)          # Finding impeller outlet velocities     # MSG: Function need to be updated with classes


            """ The following five lines are made to find porperties of the slip for a given blade number and blade angle.
                    None of the resulting variables are applied in further calculations, but could easily be. They are included
                    for further development and demonstration of slip.  """
            Ctheta2ideal = Compressor.U2 - Ctheta2m * np.tan((np.abs((IterationMatrix.beta2BArr[ib])[0])))   
            CTheta2Real = sigma * Ctheta2ideal         
            Cslip1 = Ctheta2ideal - CTheta2Real
            beta2flow = np.rad2deg(np.arctan((Cslip1 + Cm2m * np.tan(np.abs(IterationMatrix.beta2BArr[ib]))) / Cm2m))
            dh0SlipCorrected =  Compressor.U2 * (CTheta2Real) - Compressor.U1t * Compressor.Ctheta1   # Alternative for finding fluid enthalpy change, for comparison




            """ The while loop iterates for a given combination of blade number and blade angle. It iterates the efficiency demand
                    while the estimated pressure ratio is not within the given tolerance of the desired pressure ratio. It is updated in 
                    the function called pressureOverUnderEstimate. The guess efficiency is increased or decreased varying if the 
                    oressure ratio is under- or overestimated. """
            while ( (trueFalseCheck == False) and (etaStage < Compressor.etaUpperLimit and etaStage > Compressor.etaLowerLimit) ):
                
                # -------Fluid------------ Updating work with the new isentropic efficiency -------------------
                Wx = Compressor.dh0s / etaStage                         # Specific work [J/kg/K]
                workError = np.abs(Wx - dh0SlipCorrected) / Wx          # Comparing the two methods of finding enthalpy change
                # Wx = dh0SlipCorrected                                 # Work set to be given by slip corrected euler equation 
                

                """ Impeller outlet calculation. Stagnation properties denoted by zero. Isentropic relations applied on next four rows. """
                T02m = InletConditions.T00 + Wx * (Fluid.k - 1) / (Fluid.k * Fluid.R)                                                     # Stagnation exit temperature [K]     , from dh=cp*Dt   
                P02m = InletConditions.P00 * ((Compressor.dh0s * (Fluid.k - 1)  / (Fluid.k * Fluid.R * InletConditions.T00)) + 1) ** (Fluid.k / (Fluid.k - 1))     # Exit Stagnation Pressure, isentropic [Pa]      
                T2m = T02m - (Fluid.k - 1) / (2 * Fluid.k * Fluid.R) * (C2 ** 2)             # Exit temperature [K]                            from stagnation temperature
                P2m = P02m / ((T02m / T2m) ** (Fluid.k / (Fluid.k - 1)))                                # Exit pressure [Pa]                              from stagnation pressure, isentropic relation
                rho2m = P2m / (T2m * Fluid.R)                                                                       # Exit density [kg/m^3]                           from ideal gas law
                A2 = InletConditions.mdot / (rho2m * Cm2m)                                                                    # Exit area [m^2]                                 from continuity
                b2 = A2 / (math.pi * Compressor.D2)                                                                                        # Impeller exit cylinder height [m]               from geometry
                M2 = Compressor.U2 / np.sqrt(Fluid.k*Fluid.R*T02m)                                                   # Impeller exit blade mach number                 from definition
                

                # ------------------- Finding diffuser properties -------------------
                P3, P03, C3 = geometry_system_functions.diffuserFlow(P2 = P2m, P02 = P02m, rho2 = rho2m, C2 = C2, CpD = Compressor.CpD, AR = Compressor.AR)   # Finding diffuser properties


                # ------------------- Overall performance -------------------
                etaiterate, Prest, PressureTestOuterLoop = geometry_system_functions.systemTotalPerformance(P03 = P03, T02 = T02m, U2 = Compressor.U2, T1 = Compressor.T1, workInputCoeff = workInputCoeff, P00 = InletConditions.P00, T00 = InletConditions.T00, k = Fluid.k, Cp = Fluid.Cp, Pr = Compressor.Pr)
                
                """ Checking if iteration conditions are sustained. the iterated efficiency is only accepted if
                        the pressure estimate error is less than the iteration tolerance and if the iterated
                        efficiency is within the reasonable limits etaUpperLimit and etaLowerLimit """
                if (etaiterate > Compressor.etaUpperLimit or etaiterate < Compressor.etaLowerLimit) or ( etaStage > Compressor.etaUpperLimit or etaStage < Compressor.etaLowerLimit ):
                    trueFalseCheck = False
                elif abs(PressureTestOuterLoop) > Compressor.iterTol:
                    checkOuterLoop = False
                    trueFalseCheck = False
                else:
                    trueFalseCheck = True

                # trueFalseCheck = True


                if etaStage <= Compressor.etaLowerLimit or etaStage >= Compressor.etaUpperLimit:
                    break
                
                """ Updating values from nan's if iteration criterion is upheld"""
                if  ( trueFalseCheck == True ) and Compressor.U2 < Compressor.bladeVelUpperLimit:
                    IterationMatrix.etaMat[ib, iz] = etaiterate
                    IterationMatrix.pressErrorMat[ib, iz] = abs(PressureTestOuterLoop)
                    IterationMatrix.MachExitMat[ib, iz] = M2
                    IterationMatrix.b2Mat[ib, iz] = b2
                    IterationMatrix.PrestMat[ib, iz] = Prest
                    IterationMatrix.WxMat[ib, iz] = Wx
                    IterationMatrix.dh0SlipCorrMAt[ib, iz] = dh0SlipCorrected
                    IterationMatrix.c2Mat[ib, iz] = C2
                    IterationMatrix.c2mMat[ib, iz] = Cm2m
                    IterationMatrix.Ctheta2Mat[ib, iz] = Ctheta2m
                    IterationMatrix.VslipMat[ib, iz] = Cslip1
                    IterationMatrix.sigmaMat[ib, iz] = sigma
                    IterationMatrix.beta2flowMat[ib, iz] = beta2flow



                """ Updating efficiency for next iteration"""
                etaStage = pressureOverUnderEstimate(PressureTestOuterLoop, etaStage)   # MSG: Function needs to be updated with classes
        
        
            #     debug = 1       # Breakpoint for debugging
            # debug =1            # Breakpoint for debugging
                


def print_geometry_things(Compressor, IterationMatrix):
    # ------------- Preparing data for nice plotting ------------ 
    """ Plotting functions plotCompressor.py, plotFlow.py, plotSystem.py, plotVelocities.py and plotText.py are 
            called from logic.py but prepared here since data is produced here."""


    """ making nice text to go with plots """
    maxEta = np.nanmax(IterationMatrix.etaMat)
    minEta = np.nanmin(IterationMatrix.etaMat)
    countTrue = np.count_nonzero(~np.isnan(IterationMatrix.etaMat))
    totalCases = len(IterationMatrix.ZBarr) * len(IterationMatrix.beta2BArr)

    text1 = "\n" \
            r"Valid cases: " + str(countTrue) + '/'  +str(totalCases) + "\n" \
            r"$\eta _{max}$ = " + str(round(maxEta, 3)) + "   " \
            r"$\eta _{min}$ = " + str(round(minEta, 3)) + "\n" \
            r"$W_{1t}= $" + str(round(Compressor.W1t, 3)) + "m/s   " \
            r"$C_{m1}= $" + str(round(Compressor.Cm1, 3)) + "m/s"+"\n" \
            r"$U_{1t}= $" + str(round(Compressor.U1t, 3)) + "m/s"+"   " \
            r"$M_{1f}= $" + str(round(Compressor.M1, 3)) + "\n" \
            r"$U_{2t}$ = "  +str(round(Compressor.U2, 3)) + "m/s"+"   " \
            r"$U_{2t,crit}= $" + str(round(Compressor.U2Crit, 3)) +"m/s"+ "\n" \

    text2 = "\n" \
            r"$r_{t1}$ = " +str(round(Compressor.r1, 3)) + r"$m$" +" \n" \
            r"$r_{h1}$ = " +str(round(Compressor.rh1, 3)) + r"$m$" +" \n" \
            r"$r_{t2}$ = " +str(round(Compressor.r2, 3)) + r"$m$" +" \n" \
            r"$N_{crit}$ = " +str(round(Compressor.Ncrit, 3)) + r"$rpm$" +" \n" \
            
    text3 = "\n" \
            r"Desired $PR^*$: " + str(Compressor.Pr) +  "\n" \
            r"$\frac{r_h}{r_1}{*}$ = "  +str(Compressor.rhDivr1) + r" $\frac{m}{m} $" +"    " \
            r"$\frac{r_1}{r_2}{*}$ = " + str(Compressor.r1Divr2) + r" $\frac{m}{m} $" +"\n" \
            r"$N^{*}$ = " +str(round(Compressor.Ndes, 1)) + r"$rpm$" +" \n" \
            r"Proposed efficiency $\eta *=$  " + str(Compressor.etaStage0) + "\n" \
            r"Exit swirl number $\lambda_2 *$: " + str(Compressor.lambda20) +  "\n" \
            r"Tolerance*: " + str(round(Compressor.iterTol, 3)) + "   " \

    """ Preparing input to plottting functions"""
    """ Common input for all plotting functions"""
    systemVar = [Compressor.etaStage0, Compressor.lambda2, Compressor.iterTol, IterationMatrix.ZBarr, IterationMatrix.beta2BArr]
    designParam = [Compressor.r1, Compressor.r2, Compressor.rh1, Compressor.Ncrit, Compressor.Ndes]
    flowVar = [Compressor.Pr, Compressor.W1t, Compressor.Cm1, Compressor.U1t, Compressor.U2, Compressor.U2Crit]

    """ One Z(...)-array goes to one plot function. For example: Zflow goes to plotFlow.py, Zsystem goes to plotSystem.py etc."""
    Zflow = [IterationMatrix.Ctheta2Mat, IterationMatrix.pressErrorMat, IterationMatrix.etaMat, IterationMatrix.MachExitMat, IterationMatrix.PrestMat, IterationMatrix.WxMat, IterationMatrix.dh0SlipCorrMAt]     # Goes to plotFlow.py
    Zcompressor = [IterationMatrix.b2Mat, IterationMatrix.VslipMat, IterationMatrix.beta2flowMat]                                                               # Goes to plotCompressor.py
    Zsystem = [IterationMatrix.WxMat, IterationMatrix.dh0SlipCorrMAt, IterationMatrix.sigmaMat, IterationMatrix.etaMat, IterationMatrix.PrestMat]                                                 # Goes to plotSystem.py
    Zvelocities = [IterationMatrix.c2Mat, IterationMatrix.Ctheta2Mat, IterationMatrix.c2mMat, IterationMatrix.MachExitMat, IterationMatrix.VslipMat]                                            # Goes to plotVelocities.py
    text = [text1, text2, text3]                                                                                # Goes to plotText.py

    """ Plot of relative velocities """
    # fig, ax11 = plt.subplots()
    # plt.plot(Compressor.Cm1i, W1ti, '--g')
    # ax11.tick_params(axis='both', which='major', labelsize=10)
    # plt.xlabel(r"$C_{m1}  [m/s]$",fontsize=12)
    # plt.ylabel(r'$W_{1t} [m/s]$', fontsize=12)
    # plt.title("Minimisation of W1t")
    # plt.grid()

    print('Critical RPM: ' + str(round(Compressor.Ncrit, 2)) )
    print('Applied RPM: ' + str(round(Compressor.Ndes, 2)) )
    print('Rotational velocity: ' + str(round(Compressor.U2, 2)))
    print('rh: ' + str(round(Compressor.rh1, 3)) + 'm' )
    print('r1: ' + str(round(Compressor.r1, 3)) + 'm')
    print('r2: ' + str(round(Compressor.r2, 3)) + 'm' )

    if countTrue == 0:
        raise Exception(f"combination of rh/r1, r1/r2 and N gave zero valid cases. RPM to low to achieve desired pressure ratio!")

    print('Geometry.py successfully run. \n')