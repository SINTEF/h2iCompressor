"""
This Python script is divided into two main steps for determination of the characteristics of a single stage compressor:
    1. Optimize the compressor stage dimensions at the design point for a range of blade numbers and blade angles
    2. Obtain the off-design performance for a single blade number and blade angle

Multiple fluids can be selected. However, all calculations assume ideal gas behaviour.
The necessary input properties are set in three toml files.

Author(s): Petter Resell (summer intern, 2024), Martin Spillum Grønli (SINTEF Energy Research, 2025), Åsmund Ervik (SINTEF Energy Research, 2025)
"""

import marimo
import types

from matplotlib import pyplot as plt    
# Set plot parameters
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams.update({'font.size': 12})

__generated_with = "0.17.7"
app = marimo.App()

settings = types.SimpleNamespace()
#inline settings.py
settings.Fluid = Fluid
settings.InletConditions = InletConditions
settings.Compressor = Compressor
settings.IterationMatrix = IterationMatrix

#inline geometry.py

#inline plot_optimized_geometry.py

#inline off_design_performance.py

#inline plot_off_design_performance.py


with app.setup:
    # Select working fluid, 'h2', 'air' or 'r-134a'
    fluid_name = 'air'
    # Select the compressor under study, "japikse", "lowe", "LUT", "NASA", NREC", "default"
    compressor_name = "LUT"

    # Load and make instances of classes for fluid, inlet conditions, compressor and iteration matrices
    Fluid = settings.Fluid(fluid_name = fluid_name)
    InletConditions = settings.InletConditions(fluid_instance = Fluid, compressor_name = compressor_name)
    Compressor = settings.Compressor(fluid_instance = Fluid, inlet_conditions_instance = InletConditions, compressor_name = compressor_name)
    IterationMatrix = settings.IterationMatrix(compressor_instance = Compressor)


@app.cell
def ondesign():
    inducer_optimization(Compressor.Ndes, Fluid, InletConditions, Compressor)
    impeller_optimization(Compressor, InletConditions, Fluid, IterationMatrix)
    plot_geometry(Compressor, InletConditions, IterationMatrix)
    plt.show()


@app.cell
def offdesign():
    results_off_design = off_design_performance(Compressor, Fluid, InletConditions, IterationMatrix)
    plot_off_design(results_off_design, Compressor, InletConditions)
    plt.show()


if __name__ == "__main__":
    app.run()
