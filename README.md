# Aochi2018

This repertory provides the programs calculating the depth-dependent stress accumulation on faults for earthquake simulations.

The formulations are based on the following papers, using Coulomb friction and Mohr circle. 

1. Aochi, H., Dynamic asymmetry of normal and reverse faults due to constrained depth-dependent stress accumulation, Geophys. J. Int., 215, 2134-3243, doi:10.1093/ggy407, 2018.

2. Aochi, H. & K. Tsuda, Dynamic rupture simulatios based on depth-dependent stress accumulation, submitted to Geophys. J. Int., 2022.

# Programs (Aochi, 2018)

general4_distrib.f: Stress condition for a normal fault. 

general4r_distrib.f: Stress condition for a reverse fault. 

## Compile: 

gfortran general4_distrib.f

## Default setting: 

  fs = 0.6 : static frictional coefficient
  
  fd = 0.48 : dynamic frictional coefficient
  
  cohesive = 5.0 (MPa): cohesive force
  
  dip = 45.0: fault dip
  
  ( t-value = 1.0 implicitly: defining how large Mohr circle is between two friction lines according to Aochi and Ulrich, BSSA, 2015.)

## Input: 

None (all the parameters are to set in the main program inside)

## Output: "test.out" (zdep, s1, s3, fd, tau0, tp, tr, sn0)

  zdep (km) : variable.

  s1 : maximum principal stress = vertical: calculated.

  s3 : minimum prinicipal stress = horizontal: calculated. 

  fd : dynamic frictional coefficient : parameter
  
  tau0 : initial shear stress for a given fault : calculated
  
  tp : peak strength for a given fault: calculated
  
  tr : residual strengh for a given fault: calculated
  
  sn0 : initial normal stress for a given fault: calculated


# Programs (Aochi & Tsuda, 2022)
