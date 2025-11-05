"""
This Python script computes the dimensions and performance of a multistage centrifugal compressor.

TO-DO:
        - Inlet velocity, pressure and temperature are all dependent on mass flow. This variation must be accounted 
            for when feeding into the next compressor stage. This is not done now.

Author(s): Martin Spillum Grønli (SINTEF Energy Research, 2025)
"""


# Import
import numpy as np
import matplotlib.pyplot as plt


def plot_multistage_pressure_rise(stage_compressor, stage_results):
    """ Plot pressure as a function of mass flow for each stage """
    
    plt.figure()
    P_inlet = 25.1e5        # Inlet stagnation pressure of first stage [Pa]
    for i in range(len(stage_results)):
        for result in stage_results[i]:
            if i == 0:
                plt.plot(result['mdoto'], np.ones(len(result['mdoto'])) * P_inlet / 1e5, label = 'Inlet')
            plt.plot(result['mdoto'], result['P3o'] / 1e5, label = 'Stage ' + str(i + 1))
            break           # Plot pressure ratio only for first rotational speed
    
    plt.xlabel(r'$\dot{m}$' + ' ' + '[kg/s]', fontsize = 12)
    plt.ylabel(r'$P_\mathrm{outlet}$ [bar]', fontsize = 12)
    plt.legend()


def main():
    """ Main function that finds and plots compressor geometry and compressor off-design performance """
    
    fluid_name = 'h2'                   # Select working fluid, 'h2', 'air' or 'r-134a'
    
    # Select the compressor under study, "japikse", "lowe", "LUT", "NASA", NREC", "default"
    compressor_name = "NREC"
    
    # Import modules
    from h2iCompressor import settings
    from h2iCompressor import geometry
    from h2iCompressor import plot_optimized_geometry
    from h2iCompressor import off_design_performance
    from h2iCompressor import plot_off_design_performance
    
    from matplotlib import pyplot as plt    

    # Initialize fluid
    Fluid = settings.Fluid(fluid_name = fluid_name)
    
    n_stages = 3
    stage_compressor = []
    stage_inlet_conditions = []
    stage_iteration_matrices = []

    for i in range(n_stages):
        stage_inlet_conditions.append(settings.InletConditions(fluid_instance = Fluid, compressor_name = compressor_name))
        stage_compressor.append(settings.Compressor(fluid_instance = Fluid, inlet_conditions_instance = stage_inlet_conditions[i], compressor_name = compressor_name))
        stage_iteration_matrices.append(settings.IterationMatrix(compressor_instance = stage_compressor[i]))

        if i != 0:
            stage_inlet_conditions[i].P00 = P_outlet
            stage_inlet_conditions[i].T00 = T_outlet

        # Run geometry calculations
        print('\nCalculating geometry...') 
        geometry.inducer_optimization(stage_compressor[i].Ndes, Fluid, stage_inlet_conditions[i], stage_compressor[i])
        geometry.impeller_optimization(stage_compressor[i], stage_inlet_conditions[i], Fluid, stage_iteration_matrices[i])
        plot_optimized_geometry.plot_geometry(stage_compressor[i], stage_inlet_conditions[i], stage_iteration_matrices[i])
        
        P_outlet = stage_inlet_conditions[i].P00 * stage_compressor[i].Pr       # MSG: Change the last factor to the actual pressure ratio predicted by off-design script? Also this is dependent on mass flow
        T_outlet = stage_inlet_conditions[i].T00                                # MSG: Assumes perfect intercooling


        #plt.show()
    
    
    # Compute off-design performance of each stage
    stage_results = []
    for i in range(n_stages):
        results_off_design = off_design_performance.off_design_performance(stage_compressor[i], Fluid, stage_inlet_conditions[i], stage_iteration_matrices[i])
        stage_results.append(results_off_design)
        plot_off_design_performance.plot_off_design(results_off_design, stage_compressor[i], stage_inlet_conditions[i])
    plot_multistage_pressure_rise(stage_compressor, stage_results)
    
    plt.show()
    

if __name__ == '__main__':
    main()
