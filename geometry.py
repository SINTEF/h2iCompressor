"""
The following Python code takes a set of variables as inputs and calculates the geometry of an inducer, 
impeller and diffuser for a centrifugal compressor. This is done for a range of blade numbers and blade angles.

The method follows that in Japikse (1996).

References: Japikse, David (1996), Centrifugal Compressor Design and Performance, page. 6-4
            Galvas, Michael R. (1973). Fortran program for predicting off-design performance of centrifugal compressors, https://ntrs.nasa.gov/citations/19740001912 
            Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244

TO-DO:
        - Make it possible to set a maximum efficiency and then compute the maximal pressure rise for that efficiency and rotational speed
        - Rotor efficiency is now set constant. Is there a more precise way to calculate it?

Author(s): Petter Resell (summer intern, 2024), Martin Spillum Grønli (SINTEF Energy Research, 2025)
"""


# Import
import numpy as np


def inducer_optimization(N, Fluid, InletConditions, Compressor):
    """ Optimize inducer by minimizing relative velocity at inducer exit """
    
    # Make arrays for all possible values of inlet speed Cm1
    Ctheta1i = InletConditions.Cm1i * np.tan(np.radians(InletConditions.alpha1))                # Inducer absolute tangential velocity [degrees]        Trigonometry
    C1i = (Ctheta1i ** 2 + InletConditions.Cm1i ** 2) ** 0.5                                    # Inducer absolute velocity C1 [m/s]                    Pythagoras theorem from velocity triangle 
    T1i = InletConditions.T00i - (C1i ** 2) / (2 * Fluid.Cp)                                    # Inducer temperature [K]                               Stagnation temperature relation             
    
    M1i = C1i / ((Fluid.k * Fluid.R * T1i) ** 0.5)                                              # Inducer Mach number [-]                               By definition
    P1i = InletConditions.P00 * (T1i / InletConditions.T00) ** (Fluid.k / (Fluid.k - 1))        # Inducer pressure [Pa]                                 Isentropic relation                        
    rho1i = P1i / (Fluid.R * T1i)                                                               # Inducer density  [kg/m^3]                             Assuming ideal gas                          
    A1i = InletConditions.mdot / (rho1i * InletConditions.Cm1i * (1 - InletConditions.B1))      # Inducer flow area [m^2]                               Continuity                                                        
    if Compressor.optimize_rh_and_r1:        # If rh/r1 is given
        rt1i = (A1i / (np.pi * (1 - Compressor.rhDivr1 ** 2))) ** 0.5                           # Inducer tip radius [m]                                Geometry
    else:   # If rh is given and not rh/r1
        rt1i = (A1i / np.pi + Compressor.rh ** 2) ** 0.5
    U1ti = 2 * np.pi * rt1i * N / 60                                                            # Inducer blade tip speed [m/s]
    W1ti = (InletConditions.Cm1i ** 2 + (U1ti - Ctheta1i) ** 2) ** 0.5                          # Inducer relative velocity [m/s]
    
    Compressor.W1ti = W1ti                                                                      # Store relative velocities for plotting

    # Get the index of the minimum relative velocity
    w1ti_index_min = np.argmin(W1ti)                                                            # Index of smallest relative velocity [m/s]
    
    # Set inducer properties to the values that minimize the relative velocity
    Compressor.r1 = rt1i[w1ti_index_min]                                                        # Inducer tip radius [m]
    if Compressor.optimize_rh_and_r1:
        Compressor.rh = Compressor.rhDivr1 * Compressor.r1
    Compressor.Ctheta1 = Ctheta1i[w1ti_index_min]                                               # Inducer angular velocity [m/s]
    Compressor.C1 = C1i[w1ti_index_min]                                                         # Inducer velocity [m/s]
    Compressor.T1 = T1i[w1ti_index_min]                                                         # Inducer temperature [K]
    Compressor.M1 = M1i[w1ti_index_min]                                                         # Inducer Mach number [-]
    Compressor.P1 = P1i[w1ti_index_min]                                                         # Inducer pressure [Pa]
    Compressor.rho1 = rho1i[w1ti_index_min]                                                     # Inducer density [kg/m3]
    Compressor.A1 = A1i[w1ti_index_min]                                                         # Inducer flow area [m2]
    Compressor.U1t = U1ti[w1ti_index_min]                                                       # Inducer tip speed [m/s]
    Compressor.Cm1 = InletConditions.Cm1i[w1ti_index_min]                                       # Inducer meridional velocity [m/s]
    Compressor.beta1 = np.degrees(np.arctan((Compressor.U1t - Compressor.Ctheta1) / Compressor.Cm1))    # Inducer relative velocity angle [deg]
    Compressor.W1t = W1ti[w1ti_index_min]                                                       # Inducer relative velocity [m/s]
    Compressor.omega = Compressor.U1t / Compressor.r1                                           # Angular velocity [rad/s]   


def impeller_optimization(Compressor, InletConditions, Fluid, IterationMatrix):
    """ For all combinations of blade number and blade angle, optimize impeller dimensions and obtain pressure ratio for a given rotor efficiency, diffuser area ratio and diffuser efficiency """    

    for iz in range(len(IterationMatrix.ZB)):
        for ib in range(len(IterationMatrix.beta2B)):
            blade_number = IterationMatrix.ZB[iz]
            blade_angle = np.rad2deg(IterationMatrix.beta2B[ib])
            print('\nBlade number =', blade_number, 'Blade angle =', blade_angle)
            
            sigma = IterationMatrix.sigmaWiesner[ib, iz]                         # Slip factor for current blade number and blade angle
            etaStage = 0.9      # Initial guess for stage efficiency 
            Pr_stage = 0        # Initial value of stage pressure ratio
            eps = 0.0001        
            while abs(Pr_stage / Compressor.Pr - 1) > eps:
                dh0s = ((Fluid.k * Fluid.R * InletConditions.T00) / (Fluid.k - 1)) * (Compressor.Pr ** ((Fluid.k - 1) / Fluid.k) - 1)       
                Wx = dh0s / etaStage                                                # Specific work [J/kg/K]
                
                T02m = InletConditions.T00 + Wx * (Fluid.k - 1) / (Fluid.k * Fluid.R)        # Exit stagnation temperature [K]          from dh=cp*Dt   
                mu = sigma * Compressor.lambda2 / (Compressor.lambda2 - np.tan((IterationMatrix.beta2B[ib])))
                U2 = ((Compressor.U1t * Compressor.Ctheta1 + Wx) / mu) ** (1 / 2)   # Impeller tip velocity [m/s]                       from work input coefficient
                D2 = 60 * U2 / (np.pi * Compressor.Ndes)
                r2 = D2 / 2
                Ctheta2m = mu * U2
                Cm2m = Ctheta2m / Compressor.lambda2
                C2 = (Ctheta2m ** 2 + Cm2m ** 2) ** 0.5 
                T2m = T02m - (Fluid.k - 1) / (2 * Fluid.k * Fluid.R) * C2 ** 2      # Exit static temperature [K]              
                M2 = U2 / np.sqrt(Fluid.k * Fluid.R * T2m)                          # Impeller exit blade mach number              
                P02m = InletConditions.P00 * ((Wx * Compressor.eta_rotor * (Fluid.k - 1)  / (Fluid.k * Fluid.R * InletConditions.T00)) + 1) ** (Fluid.k / (Fluid.k - 1))     # Exit Stagnation Pressure, isentropic [Pa]                   
                P2m = P02m / ((T02m / T2m) ** (Fluid.k / (Fluid.k - 1)))            # Exit pressure [Pa]                                from stagnation pressure, isentropic relation
                rho2m = P2m / (T2m * Fluid.R)                                       # Exit density [kg/m^3]                             from ideal gas law
                A2m = InletConditions.mdot / (rho2m * Cm2m)                         # Exit area [m^2]                                   from continuity
                b2 = A2m / (np.pi * D2)                                             # Impeller exit cylinder height [m]                 from geometry
                
                CpDi = 1 - 1 / Compressor.AR ** 2   
                CpD = Compressor.etad * CpDi
                
                p5 = P2m + CpD * (P02m - P2m)                                       # From Bernoulli's eq. for a vaneless diffuser

                etaStage = ((p5 / InletConditions.P00) ** ((Fluid.k - 1) / Fluid.k) - 1) / (T02m / InletConditions.T00 - 1)
                Pr_stage = p5 / InletConditions.P00

            r_inlet_rms = np.sqrt(1 / 2 * (Compressor.rh ** 2 + Compressor.r1 ** 2))
            if r_inlet_rms / r2 > np.exp(- 8.16 * np.cos(IterationMatrix.beta2B[ib]) / IterationMatrix.ZB[iz]):     # Add correction factor if this becomes and issue. Some authors use r1/r2 for LHS of condition
                raise ValueError("Wiesner slip factor correlation is not valid for the radius ratio r_inlet_rms / r2 = " + str(r_inlet_rms / r2) + ' > ' + str(np.exp(- 8.16 * np.cos(IterationMatrix.beta2B[ib]) / IterationMatrix.ZB[iz])))

            # Save conditions for current blade number and blade angle
            IterationMatrix.eta[ib, iz] = etaStage
            IterationMatrix.Pr[ib, iz] = Pr_stage
            IterationMatrix.r2[ib, iz] = r2  
            IterationMatrix.b2[ib, iz] = b2
            IterationMatrix.U2[ib, iz] = U2
            IterationMatrix.Wx[ib, iz] = Wx
            IterationMatrix.P02m[ib, iz] = P02m
            IterationMatrix.P2m[ib, iz] = P2m
            IterationMatrix.p5[ib, iz] = p5
            IterationMatrix.T02m[ib, iz] = T02m
            IterationMatrix.T2m[ib, iz] = T2m
            IterationMatrix.M2[ib, iz] = M2
            IterationMatrix.c2[ib, iz] = C2
            IterationMatrix.cm2m[ib, iz] = Cm2m
            IterationMatrix.Ctheta2[ib, iz] = Ctheta2m
            IterationMatrix.sigma[ib, iz] = sigma

        # Set design point dimensions         
        design_point_blade_angle_index = np.where(IterationMatrix.beta2B == Compressor.bladeAngle)[0][0]   
        design_point_blade_number_index = np.where(IterationMatrix.ZB == Compressor.bladeNumber)[0][0]
        Compressor.r2 = IterationMatrix.r2[design_point_blade_angle_index][design_point_blade_number_index]
        Compressor.D2 = 2 * Compressor.r2 
        Compressor.U2 = IterationMatrix.U2[design_point_blade_angle_index][design_point_blade_number_index]
        Compressor.b2 = IterationMatrix.b2[design_point_blade_angle_index][design_point_blade_number_index]
        Compressor.etaStage = IterationMatrix.eta[design_point_blade_angle_index][design_point_blade_number_index]
        Compressor.dh0s = dh0s

    print('\nGeometry successfully calculated')  

    print('\nDesign point:')
    print('\tRotational speed = ' + str(round(Compressor.Ndes, 2)) + ' rpm' )
    print('\trh = ' + str(round(Compressor.rh, 10) * 1000) + ' mm' )
    print('\tr1 = ' + str(round(Compressor.r1, 10) * 1000) + ' mm')
    print('\tr2 = ' + str(round(Compressor.r2 , 10) * 1000) + ' mm' )
    print('\tb2 = ' + str(round(Compressor.b2 , 10) * 1000) + ' mm' )  
    print('\tetaStage = ' + str(round(Compressor.etaStage , 10)))  


if __name__ == '__main__':
    pass