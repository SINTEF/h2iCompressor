"""
This Python script is divided into two main steps for determination of the characteristics of a single stage compressor:
    1. Optimize the compressor stage dimensions at the design point for a range of blade numbers and blade angles
    2. Obtain the off-design performance for a single blade number and blade angle

Multiple fluids can be selected. However, all calculations assume ideal gas behaviour.
The necessary input properties are set in three toml files.

Author(s): Petter Resell (summer intern, 2024), Martin Spillum Grønli (SINTEF Energy Research, 2025)
"""


def main():
    """ Main function that finds and plots compressor stage geometry and compressor stage off-design performance """
    
    # Select working fluid, 'h2', 'air' or 'r-134a'
    fluid_name = 'air'              
    
    # Toml files with necessary input properties
    path_to_fluid_properties_toml = './properties/fluid_' + fluid_name + '.toml'
    path_to_inlet_toml = './properties/inlet_conditions_lut.toml'
    path_to_compressor_toml = './properties/compressor_lut.toml'
   
    # Import 
    import settings
    import geometry
    import plot_optimized_geometry
    import off_design_performance
    import plot_off_design_performance
    
    from matplotlib import pyplot as plt    
    
    # Set plot parameters
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'font.size': 12})

    # Load and make instances of classes for fluid, inlet conditions, compressor and iteration matrices
    Fluid = settings.Fluid(path_to_fluid_properties_toml = path_to_fluid_properties_toml)
    InletConditions = settings.InletConditions(fluid_instance = Fluid, path_to_inlet_toml = path_to_inlet_toml)
    Compressor = settings.Compressor(fluid_instance = Fluid, inlet_conditions_instance = InletConditions, path_to_compressor_toml = path_to_compressor_toml)
    IterationMatrix = settings.IterationMatrix(compressor_instance = Compressor)
    
    print('\nCalculating geometry...') 
    geometry.inducer_optimization(Compressor.Ndes, Fluid, InletConditions, Compressor)
    geometry.impeller_optimization(Compressor, InletConditions, Fluid, IterationMatrix)
    plot_optimized_geometry.plot_geometry(Compressor, InletConditions, IterationMatrix)
    plt.show()

    print('\nCalculating off-design performance...')
    results_off_design = off_design_performance.off_design_performance(Compressor, Fluid, InletConditions, IterationMatrix)
    plot_off_design_performance.plot_off_design(results_off_design, Compressor, InletConditions)
    plt.show()
       

if __name__ == '__main__':
    main()