"""
This Python script is divided into two main steps for determination of the characteristics of a single stage compressor:
    1. Optimize the compressor stage dimensions at the design point for a range of blade numbers and blade angles
    2. Obtain the off-design performance for a single blade number and blade angle

Multiple fluids and compressors can be selected. All calculations assume ideal gas behaviour.

Author(s): Petter Resell (summer intern, 2024), Martin Spillum Grønli (SINTEF Energy Research, 2025), Åsmund Ervik (SINTEF Energy Research, 2025)
"""

import marimo

__generated_with = "0.17.7"
app = marimo.App()


@app.cell
def select():
    import marimo as mo

    # Select working fluid, 'h2', 'air' or 'r-134a'
    fluid_name = mo.ui.dropdown(options=["air", "h2", "R134a"], value="air", label="Working fluid")

    # Select the compressor under study, 
    compressor_name = mo.ui.dropdown(options=["japikse", "lowe", "LUT", "NASA", "NREC", "default"], value="LUT", label="Compressor specification")
    (fluid_name,compressor_name)
    return compressor_name, fluid_name


@app.cell
def init(compressor_name, fluid_name):
    from matplotlib import pyplot as plt    
    # Set plot parameters
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'font.size': 12})

    #from h2iCompressor import settings
    import types
    settings = types.SimpleNamespace()
    #inline settings.py
    settings.Fluid = Fluid
    settings.InletConditions = InletConditions
    settings.Compressor = Compressor
    settings.IterationMatrix = IterationMatrix
    # Load and make instances of classes for fluid, inlet conditions, compressor and iteration matrices
    Fluid = settings.Fluid(fluid_name = fluid_name.value)
    InletConditions = settings.InletConditions(fluid_instance = Fluid, compressor_name = compressor_name.value)
    Compressor = settings.Compressor(fluid_instance = Fluid, inlet_conditions_instance = InletConditions, compressor_name = compressor_name.value)
    IterationMatrix = settings.IterationMatrix(compressor_instance = Compressor)

    print(f"Initialised fluid {Fluid.name} and compressor type {Compressor.name}")

    return Compressor, Fluid, InletConditions, IterationMatrix, plt


@app.cell
def ondesign(Compressor, Fluid, InletConditions, IterationMatrix, plt):

    #from h2iCompressor.geometry import inducer_optimization, impeller_optimization
    #from h2iCompressor.plot_optimized_geometry import plot_geometry
    #inline geometry.py
    #inline plot_optimized_geometry.py

    inducer_optimization(Compressor.Ndes, Fluid, InletConditions, Compressor)
    impeller_optimization(Compressor, InletConditions, Fluid, IterationMatrix)
    plot_geometry(Compressor, InletConditions, IterationMatrix)
    plt.show()
    return


@app.cell
def offdesign(Compressor, Fluid, InletConditions, IterationMatrix, plt):
    #from h2iCompressor.off_design_performance import off_design_performance
    #from h2iCompressor.plot_off_design_performance import plot_off_design
    #inline off_design_performance.py
    #inline plot_off_design_performance.py

    results_off_design = off_design_performance(Compressor, Fluid, InletConditions, IterationMatrix)
    plot_off_design(results_off_design, Compressor, InletConditions)
    plt.show()
    return


if __name__ == "__main__":
    app.run()
