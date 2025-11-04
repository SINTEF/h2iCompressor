# Centrifugal compressor code

This python code is aimed at predicting centrifugal compressor performance for a specified fluid and geometry.
The workflow is divided into two steps:
1. Predict centrifugal compressor geometry at design point using method proposed by Japikse (1996)
2. Predict off-design performance using method proposed by Galvas (1973) and adapted by Lowe (2016)

The code is fully implemented for a single stage centrifugal compressor, for which it has been partly validated so far.
Scripts for multistage compressors have been implemented, but needs further improvement.

### Disclaimer: 
There is some discrepancy between the pressure ratio calculated in step 1 compared to step 2.
Please see the pdf in documentation for more details.

### References:
Japikse, David (1996). _Centrifugal Compressor Design and Performance_, page. 6-4<br />
Galvas, Michael R. (1973). _Fortran program for predicting off-design performance of centrifugal compressors_. https://ntrs.nasa.gov/citations/19740001912<br />
Lowe, Mitchell (2016). _Design of a load dissipation device for a 100 kW supercritical CO2 turbine_. https://doi.org/10.14264/uql.2017.244<br />

Authors: Petter Resell, Martin Spillum Grønli, Åsmund Ervik; SINTEF Energy Research