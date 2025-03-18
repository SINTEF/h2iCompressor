"""
The following Python code iterates through a range of rotational speeds to deteremine compressor off-design performance.
A geometrical basis is found in geometry.py which is further built on in off-design_performance.py. Plotting of several variables of interest 
is done by utilizing plot_compessor.py, plot_system.py, plot_text.py, plot_velocities.py and pressure_test.py.
The script called settings.py is used for all parametrization of for instance fluid properties, inlet flow conditions, diffuser propertis etc. 

Authors: Petter Resell (summer intern, 2024), Martin Spillum Grønli (SINTEF Energy Research, 2024)
"""


def main():
    """ Main function that finds and plots compressor geometry and compressor off-design performance """
    fluid_name = 'air'                   # Select working fluid, 'h2' or 'air'
    path_to_fluid_properties_toml = './properties/' + fluid_name + '.toml'
    #path_to_inlet_toml = './properties/inlet_conditions_moen.toml'
    #path_to_compressor_toml = './properties/compressor_moen.toml'
    path_to_inlet_toml = './properties/inlet_conditions_lut.toml'
    path_to_compressor_toml = './properties/compressor_lut.toml'
    #path_to_inlet_toml = './properties/inlet_conditions_nrec_h2.toml'
    #path_to_compressor_toml = './properties/compressor_nrec_h2.toml'
   
    # Import necessary modules
    import settings
    import geometry
    import geometry_2
    import off_design_performance
    
    import matplotlib       
    matplotlib.use('tkagg')     # Use tkagg to avoid crashing when using X11-forwarding for plotting
    from matplotlib import pyplot as plt    

    # Load and make instances of classes for fluid, inlet conditions, compressor and iteration matrices
    Fluid = settings.Fluid(path_to_fluid_properties_toml = path_to_fluid_properties_toml)
    InletConditions = settings.InletConditions(fluid_instance = Fluid, path_to_inlet_toml = path_to_inlet_toml)
    Compressor = settings.Compressor(fluid_instance = Fluid, inlet_conditions_instance = InletConditions, path_to_compressor_toml = path_to_compressor_toml)
    IterationMatrix = settings.IterationMatrix(compressor_instance = Compressor)
    
    print('\nCalculating geometry...') 
    if Compressor.optimize_inducer_geometry:
        geometry.inducer_and_impeller_calculations(Fluid, InletConditions, Compressor)
        geometry.iterate_blade_number_and_blade_angle(Compressor, InletConditions, Fluid, IterationMatrix)
        geometry.print_and_plot_geometry(Compressor, InletConditions, IterationMatrix)
    else:
        geometry_2.inducer_and_impeller_calculations(Fluid, InletConditions, Compressor)
        geometry_2.iterate_blade_number_and_blade_angle(Compressor, InletConditions, Fluid, IterationMatrix)
        geometry_2.print_and_plot_geometry(Compressor, InletConditions, IterationMatrix)
    #plt.show()



    print('\nCalculating off-design performance...')
    results, constant_eff_line_mdot, constant_eff_line_pr, constant_eff_lines = off_design_performance.off_design_performance(Compressor, Fluid, InletConditions, IterationMatrix)
    off_design_performance.plot_off_design_performance(results, constant_eff_line_mdot, constant_eff_line_pr, constant_eff_lines, Compressor, InletConditions)
    plt.show()

    #off_design_performance.off_design_performance(Compressor, Fluid, InletConditions, IterationMatrix)
    

if __name__ == '__main__':
    main()