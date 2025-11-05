"""
The following Python code defines classes which load variables for the working fluid, inlet conditions and compressor from toml files. 


References: Japikse, David (1996), Centrifugal Compressor Design and Performance, page. 6-4
            
            Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244
Author(s): Petter Resell (summer intern, 2024), Martin Spillum Grønli (SINTEF Energy Research, 2024)
"""

# Import
import numpy as np
import toml

class Fluid:
    """ Class with fluid properties """
    
    def __init__(self, fluid_name):
        self.fluid_properties = toml.load("fluids.toml")[fluid_name]    # Load properties from toml file
        self.name = fluid_name
        self.Cp = self.fluid_properties['Cp']
        self.MolarMass = self.fluid_properties['MolarMass']
        self.k = self.fluid_properties['k']
        self.viscosity = self.fluid_properties['viscosity']
        self.R = 8.314 / self.MolarMass


class InletConditions:
    """ Class with inlet conditions """
    
    def __init__(self, fluid_instance, compressor_name):
        self.inlet_conditions = toml.load("compressors.toml")[compressor_name]['inlet_conditions']               # Load inlet conditions from toml file
        self.mdot = self.inlet_conditions['mdot']
        self.P00 = self.inlet_conditions['P00']
        self.T00 = self.inlet_conditions['T00']
        self.alpha1 = self.inlet_conditions['alpha1']
        self.B1 = self.inlet_conditions['B1']
        self.rho0 = self.P00 / (fluid_instance.R * self.T00)   

        # If optimizing inducer geometry by varying r1, rh and r2
        if 'Cm1i_min' in self.inlet_conditions and 'Cm1i_max' in self.inlet_conditions and 'Cm1i_step' in self.inlet_conditions:
            self.Cm1i = np.arange(self.inlet_conditions['Cm1i_min'], self.inlet_conditions['Cm1i_max'], self.inlet_conditions['Cm1i_step']) 
            self.T00i = np.full(len(self.Cm1i), self.T00)

class Compressor:
    """ Class with compressor properties """
    
    def __init__(self, fluid_instance, inlet_conditions_instance, compressor_name):
        self.name = compressor_name
        self.compressor = toml.load("compressors.toml")[compressor_name]    # Load compressor properties from toml file
        self.impeller_properties = self.compressor['impeller_properties']
        self.impellerDensity = self.impeller_properties['impellerDensity']
        self.impellerTensileStrength = self.impeller_properties['impellerTensileStrength']
        
        # Calculating critical property of a rotating disk to use as constraint for RPM, rotational velocity or radius. Disk will break at outermost point, therefore r2 and U2.
        self.U2Crit = np.sqrt(2 * self.impellerTensileStrength / self.impellerDensity)          # Applying tensile strength of disk. Titanium used.        
        
        # Inducer parameters
        self.W1ti = None    # Inducer relative velocities [m/s] found for each Cm1
        self.r1 = None      # Inducer tip radius [m]      
        self.Ctheta1 = None # Inducer angular velocity [m/s]
        self.C1 = None      # Inducer velocity [m/s]
        self.T1 = None      # Inducer temperature [K]
        self.M1 = None      # Inducer Mach number [-]
        self.P1 = None      # Inducer pressure [Pa]
        self.rho1 = None    # Inducer density [kg/m3]
        self.A1 = None      # Inducer area [m2]
        self.U1t = None     # Inducer tip speed [m/s]
        self.Cm1 = None     # Inducer meridional velocity [m/s]
        self.beta1 = None   # Inducer relative velocity angle [deg]
        self.W1t = None     # Inducer relative velocity [m/s]
        self.omega = None   # Angular velocity [rad/s]         
        self.rh1 = None     # Hub radius [m]

        # Impeller
        self.r2 = None      # Impeller tip radius [m]
        self.D2 = None      # Impeller tip diameter [m]
        self.U2 = None      # Impeller tip speed [m/s]
        self.Ncrit = None   # Critical rotational speed [rpm]

        # Impeller exit flow parameters
        self.lambda2 = self.impeller_properties['lambda2']
        self.eta_rotor = self.impeller_properties['eta_rotor']

        self.diffuser_properties = self.compressor['diffuser_properties']
        self.etad = self.diffuser_properties['etad']
        self.AR = self.diffuser_properties['AR']
        
        self.off_design_parameters_unknown = self.compressor['off_design_parameters_unknown']
        self.No = self.off_design_parameters_unknown['No']                  # Off design speed percentage that requires performance prediction [rpm]
        self.tu = self.off_design_parameters_unknown['tu']                  # Blade thickness [m]
        self.curvet1 = self.off_design_parameters_unknown['curvet1']        # Inducer inlet tip wall curvature [m^-1]
        self.curveh1 = self.off_design_parameters_unknown['curveh1']        # Inducer inlet hub wall curvature [m^-1]
        self.x = self.off_design_parameters_unknown['x']                    # Streamline angle [degrees] from axial direction
        self.kBL = self.off_design_parameters_unknown['kBL']                # Blading loss coefficient [-]
        self.kSF = self.off_design_parameters_unknown['kSF']                # Skin friction coefficient [-]
        self.Cf = self.off_design_parameters_unknown['Cf']                  # Friction coefficient [-]
        
        self.iteration_parameters = self.compressor['iteration_parameters']
        self.beta2Bmax = self.iteration_parameters['beta2Bmax']
        self.beta2Bmin = self.iteration_parameters['beta2Bmin']
        self.bladeMin = self.iteration_parameters['bladeMin']
        self.bladeMax = self.iteration_parameters['bladeMax']

        self.parameters_to_vary = self.compressor['parameters_to_vary']
        self.Pr = self.parameters_to_vary['Pr']
        self.Ndes = self.parameters_to_vary['Ndes']

        # If optimizing inducer geometry by varying r1, rh and r2
        if 'rh' in self.parameters_to_vary:
            self.rh = self.parameters_to_vary['rh']
            self.optimize_rh_and_r1 = False
        else:
            self.rhDivr1 = self.parameters_to_vary['rhDivr1']
            self.optimize_rh_and_r1 = True

    
        self.off_design_parameters = self.compressor['off_design_parameters']
        self.bladeAngle = np.deg2rad(self.off_design_parameters['bladeAngle'])
        self.bladeNumber = self.off_design_parameters['bladeNumber']
        self.N_off_design_arr = self.off_design_parameters['N_off_design_arr']
        self.V0DivVcr_off_design_arr = self.off_design_parameters['V0DivVcr_off_design_arr']
        self.constant_eff_lines = self.off_design_parameters['constant_eff_lines']


class IterationMatrix:
    """ Class with iteration matrices """
    
    def __init__(self, compressor_instance):
        self.ZB = np.arange(compressor_instance.bladeMin, compressor_instance.bladeMax + 1, 1)                          # Array with increasing blade number
        self.beta2B = np.radians(np.arange(compressor_instance.beta2Bmax, compressor_instance.beta2Bmin + 1, 1))        # Array with decreasing (absolute) blade angle [rad]
        beta2B_flipped = self.beta2B[:, np.newaxis]                                                                     # Flipping from row to column vector to make matrix on next lines

        # Making matrices for iteration
        self.rhsExpLimit = np.exp(- 8.16 * np.cos(beta2B_flipped) / self.ZB)                                                            # Matrix for epsilon_limit from Eiesner condition           
        self.sigmaWiesner = 1 - (np.sqrt(np.cos(np.radians(beta2B_flipped))) / (self.ZB ** 0.7))                                        # Matrix for Wiesner slip factor 

        self.eta = np.array([[np.nan for _ in range(np.shape(self.rhsExpLimit)[1])] for _ in range(np.shape(self.rhsExpLimit)[0])])     # Matrix for efficiency
        self.Pr = np.copy(self.eta)                                                                                                     # Matrix for pressure estimate
        self.r2 = np.copy(self.eta)                                                                                                     # Matrix for impeller tip radius
        self.b2 = np.copy(self.eta)                                                                                                     # Matrix for impeller exit cylinder height
        self.U2 = np.copy(self.eta)                                                                                                     # Matrix for impeller rotational speed        
        self.Wx = np.copy(self.eta)                                                                                                     # Matrix for compression work
        self.P02m = np.copy(self.eta)  
        self.P2m = np.copy(self.eta)  
        self.p5 = np.copy(self.eta)  
        self.T02m = np.copy(self.eta)  
        self.T2m = np.copy(self.eta)    
        self.M2 = np.copy(self.eta)                                                                                                     # Matrix for impeller mach number
        self.c2 = np.copy(self.eta)                                                                                                     # Matrix for impeller absolute discharge velocity
        self.cm2m = np.copy(self.eta)                                                                                                   # Matrix for meridional component of impeller discharge velocity
        self.Ctheta2 = np.copy(self.eta)                                                            
        self.sigma = np.copy(self.eta)                                                                                                  # Matrix for slip factor found to be valid
