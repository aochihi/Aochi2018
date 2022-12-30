# Depth depedent stress setting for dynamic rupture simulations

This repertory provides the programs calculating the depth-dependent stress accumulation on faults for earthquake simulations.

The formulations are based on the following papers, using Coulomb friction and Mohr circle. 

1. Aochi, H., Dynamic asymmetry of normal and reverse faults due to constrained depth-dependent stress accumulation, Geophys. J. Int., 215, 2134-3243, 2018. https//doi.org/10.1093/gji/ggy407

2. Aochi, H. & K. Tsuda, Dynamic rupture simulatios based on depth-dependent stress accumulation, Geophys. J. Int., published on line, 2022. https://doi.org/10.1093/gji/ggac453

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

layer5rev_distrib.f: reverse faulting (Figure 3 and Figure 10)

layer5ss_distrib.f: strike-slip faulting (Figure 4)

## Compile

gfortran layer5rev_distrib.f

## Input:

infile = "1d.txt" or "1d_teil.txt"

1st line: Total number of layers.

2nd lines - : Depth of layer top (km), Vs (m/s), medium density (kg/m3)

## Output 1: "before_constrained.dat" for layer5rev_distrib.f

Depth (km), S_3=S_v (MPa), S_1=S_H (MPa), Delta sigma (MPa), rigidity (Pa), Delta epsilon

## Output 2: "constrained.dat" for layer5rev_distrib.f

1st-2nd lines: headers

Depth (km), S_3=S_v (MPa), S_1=S_H (MPa), Delta sigma (MPa), ridigidy (Pa), Delta epsilon

## Output 1ss: "ss_before_constrained.dat" for layer5s_distrib.f

Depth (km), S_3=S_h (MPa), S_2=S_v (MPa), S_1=S_H (Mpa), Delta sigma (MPa), rigidity (Pa), Delta epsilon

## Output 2ss: "ss_constrained.dat" for layer5ss_distrib.f

1st-2nd lines: headers

Depth (km), S_3=S_h (MPa), S_2=S_v (MPa), S_1=S_H (Mpa), Delta sigma (MPa), rigidity (Pa), Delta epsilon



