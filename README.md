# Lidar_Source

## Descripton:
This repository contains code for processing and visualizing lidar data produced by various backscatter lidars.

## System requirements:
- Python 3.X
- Numpy 1.14 or later
- Runs on 64 bit Linux, 64 or 32 bit Windows.
- IDL (optional) for a few minor visualizations/analysis codes.

## Notes:
[5/25/18]
As of this date, code has been written to process raw CPL and CAMAL data into NRB.
This project only contains the code necessary to process. The actual data files 
themselves should be located on a local computer. Paths to data files are set in 
an initialization file (initializations.py), which is contained within the source code
of this project. Although most initialization file parameters will overlap between lidars,
there will be some differences. Therefore, each lidar will have a "default" initializations
file (cpl_default_initializations.py, etc) which will serve as a reference for the required
set of initialization parameters.
