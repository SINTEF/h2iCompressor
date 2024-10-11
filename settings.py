"""
The following Python code defines classes which load variables for the working fluid, inlet conditions and compressor from toml files. 


TO-DO:
    - When defining functions in other scripts, let functions take classes as arguments instead of individual variables
    - Make sure that all variables are defined in the classes (search for ex. Compressor in other scripts)

Reference: Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244
Author(s): Petter Resell (summer intern, 2024), Martin Spillum Grønli (SINTEF Energy Research, 2024)
"""


# Import
import numpy as np
import toml


class Fluid:
    """ Class with fluid properties """
    def __init__(self, path_to_fluid_properties_toml):
        self.fluid_properties = toml.load(path_to_fluid_properties_toml)['fluid_properties']    # Load properties from toml file
        self.Cp = self.fluid_properties['Cp']
        self.MolarMass = self.fluid_properties['MolarMass']
        self.k = self.fluid_properties['k']
        self.viscosity = self.fluid_properties['viscosity']
        self.R = 8.314 / self.MolarMass


class InletConditions:
    """ Class with inlet conditions """
    def __init__(self, fluid_instance, path_to_inlet_toml):
        self.inlet_conditions = toml.load(path_to_inlet_toml)['inlet_conditions']               # Load inlet conditions from toml file
        self.mdot = self.inlet_conditions['mdot']
        self.P00 = self.inlet_conditions['P00']
        self.T00 = self.inlet_conditions['T00']
        self.alpha1 = self.inlet_conditions['alpha1']
        self.B1 = self.inlet_conditions['B1']
        self.rho0 = self.P00 / (fluid_instance.R * self.T00)   

        # If optimizing inducer geometry by varying r1, rh and r2
        if 'Cm1i_min' in self.inlet_conditions and 'Cm1i_max' in self.inlet_conditions and 'Cm1i_step' in self.inlet_conditions:
            self.Cm1i = np.arange(self.inlet_conditions['Cm1i_min'], self.inlet_conditions['Cm1i_max'], self.inlet_conditions['Cm1i_step']) # MSG: Move this to be part of Compressor under inducer? C1i is part of Compressor
            self.T00i = np.full(len(self.Cm1i), self.T00)
            self.W1ti = None      # Inducer relative velocities [m/s] found for each Cm1
        

class Compressor:
    """ Class with compressor properties """
    def __init__(self, fluid_instance, inlet_conditions_instance, path_to_compressor_toml):
        self.impeller_properties = toml.load(path_to_compressor_toml)['impeller_properties']    # Load impeller properties from toml file
        self.impellerDensity = self.impeller_properties['impellerDensity']
        self.impellerTensileStrength = self.impeller_properties['impellerTensileStrength']
        
        # Calculating critical property of a rotating disk to use as constraint for RPM/rotational velocity/radiusettings. Disk will break at outermost point, therefore r2 and U2.
        self.U2Crit = np.sqrt(2 * self.impellerTensileStrength / self.impellerDensity)          # Applying tensile strength of disk. Titanium used.        
        
        # Inducer parameters
        self.r1 = None      # Inducer tip radius [m]        # MSG: Change this to r1t?
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
        self.omega = None   # Angular velocity [rad/s]          =(2*np.pi*N/60). MSG: Inducer property?
        self.rh1 = None     # Hub radius [m]

        # Impeller
        self.r2 = None      # Impeller tip radius [m]
        self.D2 = None      # Impeller tip diameter [m]
        self.U2 = None      # Impeller tip speed [m/s]
        self.Ncrit = None   # Critical rotational speed [rpm]

        # Impeller exit flow parameters
        self.lambda2 = self.impeller_properties['lambda2']
        self.etaStage = self.impeller_properties['etaStage']
        self.etaStage0 = self.etaStage

        self.diffuser_properties = toml.load(path_to_compressor_toml)['diffuser_properties']    # Load diffuser properties from toml file
        self.etad = self.diffuser_properties['etad']
        self.AR = self.diffuser_properties['AR']
        self.CpDi = 1 - 1 / (self.AR ** 2)
        self.CpD = self.etad * self.CpDi
        
        self.off_design_parameters_unknown = toml.load(path_to_compressor_toml)['off_design_parameters_unknown']    # Load off-design parameters from toml file
        self.No = self.off_design_parameters_unknown['No']                  # Off design speed percentage that requires performance prediction [rpm]
        self.tu = self.off_design_parameters_unknown['tu']                  # Blade thickness [m]
        self.curvet1 = self.off_design_parameters_unknown['curvet1']        # Inducer inlet tip wall curvature [m^-1]
        self.curveh1 = self.off_design_parameters_unknown['curveh1']        # Inducer inlet hub wall curvature [m^-1]
        self.x = self.off_design_parameters_unknown['x']                    # Streamline angle [degrees] from axial direction
        self.V0DivVcr = np.linspace(0.8, 0.9, 300)                          # Compressor inlet absolute critical velocity ratio wrt resonance [-]       # MSG: Move to toml file. This must be updated when running with different fluids and geometries. Look at impeller mach speed, effiency, pressure ratio, etc. 
        #self.V0DivVcr = np.linspace(0.1, 0.45, 50)                         # Compressor inlet absolute critical velocity ratio wrt resonance [-]
        self.kBL = self.off_design_parameters_unknown['kBL']                # Blading loss coefficient [-]
        self.kSF = self.off_design_parameters_unknown['kSF']                # Skin friction coefficient [-]
        self.Cf = self.off_design_parameters_unknown['Cf']                  # Friction coefficient [-]
        
        
        self.iteration_parameters = toml.load(path_to_compressor_toml)['iteration_parameters']      # Load iteration parameters from toml file
        self.etaLowerLimit = self.iteration_parameters['etaLowerLimit']
        self.etaUpperLimit = self.iteration_parameters['etaUpperLimit']
        self.bladeVelUpperLimit = self.iteration_parameters['bladeVelUpperLimit']
        self.bladeVelLowerLimit = self.iteration_parameters['bladeVelLowerLimit']
        self.beta2Bmax = self.iteration_parameters['beta2Bmax']
        self.beta2Bmin = self.iteration_parameters['beta2Bmin']
        self.bladeMin = self.iteration_parameters['bladeMin']
        self.bladeMax = self.iteration_parameters['bladeMax']
        self.iterTol = self.iteration_parameters['iterTol']

        self.parameters_to_vary = toml.load(path_to_compressor_toml)['parameters_to_vary']          # Load parameters to vary from toml file
        self.Pr = self.parameters_to_vary['Pr']
        self.dh0s = ((fluid_instance.k * fluid_instance.R * inlet_conditions_instance.T00) / (fluid_instance.k - 1) ) * (self.Pr ** ((fluid_instance.k - 1) / fluid_instance.k) - 1)                     
        self.N0 = self.parameters_to_vary['N0']

        # If optimizing inducer geometry by varying r1, rh and r2
        if 'r1' in self.parameters_to_vary and 'rh' in self.parameters_to_vary and 'r2' in self.parameters_to_vary:
            self.r1 = self.parameters_to_vary['r1']
            self.rh = self.parameters_to_vary['rh']
            self.r2 = self.parameters_to_vary['r2']
            self.optimize_inducer_geometry = False
        else:
            self.rhDivr1 = self.parameters_to_vary['rhDivr1']
            self.r1Divr2 = self.parameters_to_vary['r1Divr2']
            self.optimize_inducer_geometry = True

    
        self.off_design_parameters = toml.load(path_to_compressor_toml)['off_design_parameters']    # Load off-design parameters from toml file
        self.bladeAngle = np.deg2rad(self.off_design_parameters['bladeAngle'])
        self.bladeNumber = self.off_design_parameters['bladeNumber']


class IterationMatrix:
    """ Class with iteration matrices """
    def __init__(self, compressor_instance):
        self.ZBarr = np.arange(compressor_instance.bladeMin, compressor_instance.bladeMax + 1, 1)                   # Array with increasing blade number
        self.beta2BArr = np.radians(np.arange(compressor_instance.beta2Bmax, compressor_instance.beta2Bmin, 1))     # Array with decreasing (absolute) discharge angles [rad] MSG: Is this discharge angle or blade angle?
        self.beta2BArr = self.beta2BArr[:, np.newaxis]                                                              # Flipping from row to column vector to make matrix on next lines

        """ Making matrices for iteration later. Filling with nans that are only replaced if all conditions are met. """
        self.rhsExpLimitMat = np.exp(- 8.16 * np.cos(self.beta2BArr) / self.ZBarr)                                                              # Matrix for epsilon_limit from wiesner condition           
        self.etaMat = np.array([[ np.nan for _ in range(np.shape(self.rhsExpLimitMat)[1])] for _ in range(np.shape(self.rhsExpLimitMat)[0])] )  # Matrix for efficiency. Replace by fill matrix with same shape as rhsExpLimitMAt with nans 
        self.sigmaWiesnerMat = 1 - (np.sqrt(np.cos(np.radians(self.beta2BArr))) / (self.ZBarr ** 0.7))                                          # Matrix for wiesner slip factor 
        self.sigmaMat = np.copy(self.etaMat)                                                                                                    # Matrix for slip factor found to be valid
        self.b2Mat = np.copy(self.etaMat)                                                                                                       # Matrix for impeller exit cylinder height
        self.c2Mat = np.copy(self.etaMat)                                                                                                       # Matrix for impeller absolute discharge velocity
        self.c2mMat = np.copy(self.etaMat)                                                                                                      # Matrix for meridional component of impeller discharge velocity
        self.PrestMat = np.copy(self.etaMat)                                                                                                    # Matrix for pressure estimate
        self.VslipMat = np.copy(self.etaMat)                                                                                                    # Matrix for slip velocity
        self.pressErrorMat = np.copy(self.etaMat)                                                                                               # Matrix for pressure error
        self.MachExitMat = np.copy(self.etaMat)                                                                                                 # Matrix for impeller mach number
        self.WxMat = np.copy(self.etaMat)                                                                                                       # Matrix for compression work
        self.Ctheta2Mat = np.copy(self.etaMat)                                                                                                  # Matrix for angular component of discharge velocity
        self.dh0SlipCorrMAt = np.copy(self.etaMat)                                                                                              # Matrix for work found by slip corrected euler equation
        self.beta2flowMat= np.copy(self.etaMat)                                                                                                 # Matrix for discharge flow angle found by slip relations