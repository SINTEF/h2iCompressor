"""
The following Python code iterates through a range of rotational speeds to deteremine compressor off-design performance.
A geometrical basis is found in geometry.py which is further built on in off-design_performance.py. Plotting of several variables of interest 
is done by utilizing plot_compessor.py, plot_system.py, plot_text.py, plot_velocities.py and pressure_test.py.
The script called settings.py is used for all parametrization of for instance fluid properties, inlet flow conditions, diffuser propertis etc. 

Authors: Petter Resell (summer intern, 2024), Martin Spillum Grønli (SINTEF Energy Research, 2024)
"""


def main():
    """ Main function that finds and plots compressor geometry and compressor off-design performance """
    
    fluid_name = 'h2'                   # Select working fluid, 'h2', 'air' or 'r-134a'
    
    path_to_fluid_properties_toml = './properties/fluid_' + fluid_name + '.toml'
    path_to_inlet_toml = './properties/inlet_conditions_nrec_h2.toml'
    path_to_compressor_toml = './properties/compressor_nrec_h2.toml'
   
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
    results, constant_eff_line_mdot, constant_eff_line_pr, constant_eff_lines = off_design_performance.off_design_performance(Compressor, Fluid, InletConditions, IterationMatrix)
    plot_off_design_performance.plot_off_design(results, constant_eff_line_mdot, constant_eff_line_pr, constant_eff_lines, Compressor, InletConditions)
    plt.show()
       

if __name__ == '__main__':
    main()