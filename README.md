# braggpeak_sampler

Simulation and Beamtime data analysis of the Braggpeak Sampler.

## Description
Geant4 simulation of the Braggpeak Sampler for 220 MeV Protons with no target, a RW3 target (equivalent to 5cm of water) and a LN300 target (using DCM files where 1/4 of the voxels are water and 3/4 are air).
Beamtime data analysis with root. 

## Usage
For simulation execute "cmake .." in detector/build or theory/build and then "make".
Start simulation with open OGL with ./bragg_sampler in detector/build or ./bragg_theory in theory/build.
Here the vis.mac detemines the start conditions.
For a complete run the simulation append run_p.mac to the execution command.
Change target in detector/include/DetectorConstruction.hh.

For beamtime data analysis put beamtime raw data in beamtime/MIT_05_2024/*dataset*/input.
Execute "root edep_coincidence" for full analysis in beamtime folder.

## Project status
In progress
