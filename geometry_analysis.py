import numpy as np
import toml
import settings
import geometry_new as geometry 
import off_design_performance

import matplotlib.pyplot as plt


def run_analysis(path_to_fluid_properties_toml, path_to_inlet_toml, path_to_compressor_toml, rh_r1_ratio, r1_r2_ratio, impeller_blade_number, impeller_backsweep_angle, diffuser_area_ratio):
    # Load and make instances of classes for fluid, inlet conditions, compressor and iteration matrices
    Fluid = settings.Fluid(path_to_fluid_properties_toml = path_to_fluid_properties_toml)
    InletConditions = settings.InletConditions(fluid_instance = Fluid, path_to_inlet_toml = path_to_inlet_toml)
    Compressor = settings.Compressor(fluid_instance = Fluid, inlet_conditions_instance = InletConditions, path_to_compressor_toml = path_to_compressor_toml)
    
    # Update Compressor parameters with new radius ratios, blade numbers impeller backsweep angles and diffuser area ratios before initializing the iteration matrices
    Compressor.rhDivr1 = rh_r1_ratio
    Compressor.r1Divr2 = r1_r2_ratio
    Compressor.AR = diffuser_area_ratio
    Compressor.bladeNumber = impeller_blade_number
    Compressor.bladeAngle = np.deg2rad(impeller_backsweep_angle)       # Degrees must be converted to radians since this conversion of bladeAngle is done in the initialization of the compressor class, which is already initialized
    Compressor.beta2Bmax = impeller_backsweep_angle        
    Compressor.beta2Bmin = impeller_backsweep_angle 
    Compressor.bladeMin = impeller_blade_number
    Compressor.bladeMax = impeller_blade_number             
    
    IterationMatrix = settings.IterationMatrix(compressor_instance = Compressor)
    
    # Run geometry calculations
    print('\nCalculating geometry...') 
    geometry.inducer_and_impeller_calculations(Fluid, InletConditions, Compressor)
    geometry.iterate_blade_number_and_blade_angle(Compressor, InletConditions, Fluid, IterationMatrix)

    iB = int(np.where(IterationMatrix.beta2BArr == Compressor.bladeAngle)[0])                   # Index for blade angle
    iZ = int(np.where(IterationMatrix.ZBarr == Compressor.bladeNumber)[0])                      # Index for blade number
    Compressor.eta_geometry = IterationMatrix.etaMat[iB][iZ]                                    # Efficiency at design point calculated in geometry
    #geometry.print_and_plot_geometry(Compressor, InletConditions, IterationMatrix)
    
    # Calculate off-design performance
    print('\nCalculating off-design performance...')
    results_off_design, constant_eff_line_mdot, constant_eff_line_pr, constant_eff_lines = off_design_performance.off_design_performance(Compressor, Fluid, InletConditions, IterationMatrix)
    #off_design_performance.plot_off_design_performance(results, constant_eff_line_mdot, constant_eff_line_pr, constant_eff_lines, Compressor, InletConditions)
    
    return Compressor, results_off_design


def main():
    fluid_name = 'r-134a'                   # Select working fluid, 'h2', 'r-134a' or 'air'
    path_to_fluid_properties_toml = './properties/' + fluid_name + '.toml'
    path_to_inlet_toml = './properties/inlet_conditions_japikse.toml'
    path_to_compressor_toml = './properties/compressor_japikse.toml'
    
    import matplotlib       
    matplotlib.use('tkagg')     # Use tkagg to avoid crashing when using X11-forwarding for plotting
    from matplotlib import pyplot as plt  

    # Define range of impeller blade numbers to analyze (positive integers)
    impeller_blade_numbers = np.arange(10, 30, 4)

    # Define range of impeller back sweep angles to analyze (negative integers)
    impeller_backsweep_angles = np.arange(- 40, - 39, 1)

    # Define range of radius ratios to analyze
    rh_r1_ratios = np.linspace(0.25, 0.5, 2)
    r1_r2_ratios = np.linspace(0.45, 0.50, 2)

    # Define range of diffuser area ratios to analyze
    diffuser_area_ratios = np.linspace(1.69, 3, 2)

    results = []

    print('rh/r1:', rh_r1_ratios)
    print('r1/r2:', r1_r2_ratios)
    print('Impeller backsweep angle:', impeller_backsweep_angles)
    print('Number of impeller blades:', impeller_blade_numbers)
    print('Diffuser area ratio:', diffuser_area_ratios)

    for rh_r1 in rh_r1_ratios:
        for r1_r2 in r1_r2_ratios:
            for impeller_blade_number in impeller_blade_numbers:
                for impeller_backsweep_angle in impeller_backsweep_angles:
                    for diffuser_area_ratio in diffuser_area_ratios:
                        Compressor, results_off_design = run_analysis(path_to_fluid_properties_toml, path_to_inlet_toml, path_to_compressor_toml, rh_r1, r1_r2, impeller_blade_number, impeller_backsweep_angle, diffuser_area_ratio)
                        results.append({
                            'rh_r1_ratio': rh_r1,
                            'r1_r2_ratio': r1_r2,
                            'impeller_blade_number': impeller_blade_number,
                            'impeller_backsweep_angle': impeller_backsweep_angle,
                            'diffuser_area_ratio': diffuser_area_ratio,
                            'Compressor': Compressor,
                            'results_off_design': results_off_design,
                        })
    
    # Process and display results
    process_results(results)

def process_results(results):
    # Extract data from results and plot
    rh_r1_ratios = []
    r1_r2_ratios = []
    impeller_blade_numbers = []
    impeller_backsweep_angles = []
    diffuser_area_ratios = []
    max_efficiencies = []
    max_pressure_ratios = []
    eta_geometry = []

    for result in results:
        rh_r1_ratios.append(result['rh_r1_ratio'])
        r1_r2_ratios.append(result['r1_r2_ratio'])
        impeller_blade_numbers.append(result['impeller_blade_number'])
        impeller_backsweep_angles.append(result['impeller_backsweep_angle'])
        diffuser_area_ratios.append(result['diffuser_area_ratio'])
        eta_geometry.append(result['Compressor'].eta_geometry)
        
        # Extract maximum efficiency and pressure ratio from the off-design results. MSG: This does not happen at the same mass flow rate...
        max_efficiency = max(max(result['etao']) for result in result['results_off_design'])        
        max_pressure_ratio = max(max(result['Pro']) for result in result['results_off_design'])     # Need the max pressure ratio at the design point from the off-design performance code
        
        max_efficiencies.append(max_efficiency)
        max_pressure_ratios.append(max_pressure_ratio)

    # Print a summary table
    print("\nSummary Table:")
    print("Nr. | rh/r1 ratio | r1/r2 ratio | Nr. of blades | Backsweep angle | Diffuser area ratio | Max Efficiency | Pressure Ratio | Eta geometry")
    print("-" * 115)
    design_nrs = np.arange(1, len(results) + 1)
    for design_nr, rh_r1, r1_r2, impeller_blade_number, impeller_backsweep_angle, diffuser_area_ratio, eff, pr, eta_geometry in zip(design_nrs, rh_r1_ratios, r1_r2_ratios, impeller_blade_numbers, impeller_backsweep_angles, diffuser_area_ratios, max_efficiencies, max_pressure_ratios, eta_geometry):
        print(f"{design_nr}   | {rh_r1:.2f}        | {r1_r2:.2f}        | {impeller_blade_number:.2f}         | {impeller_backsweep_angle:.2f}          | {diffuser_area_ratio:.2f}                | {eff:.4f}         | {pr:.4f} | {eta_geometry:.4f}")

    plt.figure()
    for result in results:
        Nplot = f"{result['results_off_design'][0]['N'] * 1e-3:.0f} krpm"
        plt.plot(result['results_off_design'][0]['mdoto'], result['results_off_design'][0]['Pro'], label = r'Design nr. ' + str(design_nrs[results.index(result)]))
    plt.ylabel(r'$\frac{P_{03}}{P_{00}}$', rotation = 45, fontsize = 12)
    plt.xlabel(r'$\dot{m}$' + ' ' + '[kg/s]', fontsize=12)      
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout(rect=[0, 0, 0.98, 1])
    plt.title(r'Pressure ratio at $N = $ ' + Nplot)
    
    plt.figure()
    for result in results:
        Nplot = f"{result['results_off_design'][0]['N'] * 1e-3:.0f} krpm"
        plt.plot(result['results_off_design'][0]['mdoto'], result['results_off_design'][0]['etao'], label = r'Design nr. ' + str(design_nrs[results.index(result)]))
    plt.ylabel(r'$\eta$', rotation=45, fontsize=12)
    plt.xlabel(r'$\dot{m}$' + ' ' + '[kg/s]', fontsize=12)      
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout(rect=[0, 0, 0.98, 1])
    plt.title(r'Efficiency at $N = $ ' + Nplot)

    
    # Plot correlation between max efficiency and max pressure ratio
    fig = plt.figure()
    sorted_indices = np.argsort(max_efficiencies)
    sorted_max_pressure_ratios = np.array(max_pressure_ratios)[sorted_indices]
    sorted_max_efficiencies = np.array(max_efficiencies)[sorted_indices]
    plt.plot(sorted_max_efficiencies, sorted_max_pressure_ratios)
    plt.xlabel('Max Pressure Ratio')
    plt.ylabel('Max Efficiency')
    plt.title('Correlation between Max Efficiency and Max Pressure Ratio')

    # Create a 3D scatter plot
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    scatter = ax.scatter(rh_r1_ratios, r1_r2_ratios, max_efficiencies, c=max_pressure_ratios, cmap='viridis')

    ax.set_xlabel('rh/r1 ratio')
    ax.set_ylabel('r1/r2 ratio')
    ax.set_zlabel('Max Efficiency')

    # Add a color bar
    cbar = fig.colorbar(scatter)
    cbar.set_label('Pressure Ratio')

    plt.title('Compressor Performance for Different Radius Ratios')
    plt.tight_layout()

    # Create a 2D contour plot for efficiency
    fig, ax = plt.subplots(figsize=(10, 8))
    
    unique_rh_r1 = sorted(set(rh_r1_ratios))
    unique_r1_r2 = sorted(set(r1_r2_ratios))
    
    efficiency_grid = np.zeros((len(unique_rh_r1), len(unique_r1_r2)))
    pressure_ratio_grid = np.zeros((len(unique_rh_r1), len(unique_r1_r2)))

    for i, rh_r1 in enumerate(unique_rh_r1):
        for j, r1_r2 in enumerate(unique_r1_r2):
            indices = [k for k, (rh, r1) in enumerate(zip(rh_r1_ratios, r1_r2_ratios)) if rh == rh_r1 and r1 == r1_r2]
            if indices:
                efficiency_grid[i, j] = max_efficiencies[indices[0]]
                pressure_ratio_grid[i, j] = max_pressure_ratios[indices[0]]

    contour = ax.contourf(unique_r1_r2, unique_rh_r1, efficiency_grid, cmap='viridis', levels=20)
    ax.set_xlabel('r1/r2 ratio')
    ax.set_ylabel('rh/r1 ratio')
    cbar = fig.colorbar(contour)
    cbar.set_label('Max Efficiency')

    # Add pressure ratio contour lines
    pressure_contour = ax.contour(unique_r1_r2, unique_rh_r1, pressure_ratio_grid, colors='red', levels=5)
    ax.clabel(pressure_contour, inline=True, fontsize=8, fmt='%.2f')

    plt.title('Efficiency and Pressure Ratio Contour for Different Radius Ratios')
    plt.tight_layout()
    plt.show()

    pass

if __name__ == "__main__":
    main()