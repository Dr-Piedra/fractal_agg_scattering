# fractal_agg_scattering
Generates a fractal aggregate and computes its light scattering mueller matrix


This code uses MSTM 
Mackowski, D. "MSTM Version 3.0: April 2013." (2013).

MSTM code has been set in the Repo, but it needs to be re-compiled 

to compile with gfortran:

gfortran -o mstm.exe mpidefs-serial.f90 mstm-intrinsics-v3.0.f90 mstm-modules-v3.0.f90 mstm-main-v3.0.f90

The code also uses fracagpos.f90 to generate random coordinates of monomers
This fracagpos.f90 file needs to be compiled as well and change its name to coord_generator.out

This worked for me:

gfortran fracagpos.f90
mv a.out coord_generator.out

Finally, you need to have pyquaternion module installed for python to execute the rotations; To install pyquaternion:

pip install pyquaternion

The program generate_fractal.py executes coord_generator.out and generates a file coords.data which contains the coordinates of the fractal monomers. Then, these coordinates are rotated by the Euler angles alpha,beta and gamma. Finally, the function gen_frac generates executes mstm.out using a file which includes all the angles that you want to calculate. By default the file is "scat_angles.csv". A file is created "mueller_log" in the save location (which by default is "test"). The file mueller_log and log contains the scattering phase function. 




