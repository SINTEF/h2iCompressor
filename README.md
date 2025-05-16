# Centrifugal compressor code

Python code divided into two steps:
	_1. Prediction of centrifugal compressor geometry at design point using method proposed by Japikse (1996)
	_2. Prediction of off-design performance using method proposed by Galvas (1973) and adapted by Lowe (2016)

The code is fully implemented for a single stage centrifugal compressor, for which it has been partly validated. 
Scripts for multistage compressors have been implemented, but needs further improvement.

Disclaimer: There is some discrepancy between the pressure ratio calculated in step 1 compared to step 2. Please refer to memo for more information.

References:
	Japikse, David (1996), Centrifugal Compressor Design and Performance, page. 6-4
	Galvas, Michael R. (1973). Fortran program for predicting off-design performance of centrifugal compressors, https://ntrs.nasa.gov/citations/		19740001912 
     	Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244

Authors: Petter Resell (summer intern), Martin Spillum Grønli (SINTEF Energy Research, 2024)
