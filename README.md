# Centrifugal compressor code

This python code is aimed at predicting centrifugal compressor performance for a specified fluid and geometry.
The workflow is divided into two steps:
1. Predict centrifugal compressor geometry at design point using method proposed by Japikse (1996)
2. Predict off-design performance using method proposed by Galvas (1973) and adapted by Lowe (2016)

The code is fully implemented for a single stage centrifugal compressor, for which it has been partly validated.
Scripts for multistage compressors have been implemented, but needs further work.

A live demo is available at https://sintef.github.io/h2iCompressor/

## Installation
This project uses `uv` for Python packaging, so you can get started quickly:
* If you don't already have it, you need to [install uv](https://docs.astral.sh/uv/getting-started/installation/)
* Clone this repository and enter the folder
* Run `uv sync` to install all dependencies in a virtual environment in the folder
* Use `uv run python main.py` to run the code for a single stage compressor.
* Change the values for `fluid_name` and `compressor_name` in `main.py` to run different cases
  - Not all combinations are physically sound, and you can get negative temperatures or other errors
  - Have a look inside `compressors.toml` to see the references for different compressor geometries, and the relevant working fluid(s)
* If you prefer a notebook workflow, there is a Marimo notebook corresponding to `main.py` in `notebook.py`
  - You can run this with `uv run marimo edit notebook.py`

### Disclaimers: 
There is some discrepancy between the pressure ratio calculated in step 1 compared to step 2.
Please see the pdf documentation for more details.

### References:
Japikse, David (1996). _Centrifugal Compressor Design and Performance_, page. 6-4<br />
Galvas, Michael R. (1973). _Fortran program for predicting off-design performance of centrifugal compressors_. https://ntrs.nasa.gov/citations/19740001912<br />
Lowe, Mitchell (2016). _Design of a load dissipation device for a 100 kW supercritical CO2 turbine_. https://doi.org/10.14264/uql.2017.244<br />

Authors: Petter Resell, Martin Spillum Grønli, Åsmund Ervik; SINTEF Energy Research
