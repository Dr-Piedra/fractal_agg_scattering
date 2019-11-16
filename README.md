# fractal_agg_scattering

This code generates a fractal aggregate of size parameter x, rotates the fractal aggregate at random, and computes its light scattering mueller matrix.

returns: 
- The mueller matrix
- The coordinates of the monomers


This code uses MSTM 
- http://www.eng.auburn.edu/~dmckwski/scatcodes/
- Mackowski, D. "MSTM Version 3.0: April 2013." (2013).

MSTM code has been provided the Repo, but it needs to be re-compiled 

To compile with gfortran:

gfortran -o mstm.exe mpidefs-serial.f90 mstm-intrinsics-v3.0.f90 mstm-modules-v3.0.f90 mstm-main-v3.0.f90

The code also uses fracagpos.f90 to generate random coordinates of monomers. This fracagpos.f90 file needs to be compiled as well and change its name to "coord_generator.out"

These commands worked for me:

- gfortran fracagpos.f90
- mv a.out coord_generator.out

Finally, you need to have pyquaternion module installed for python to execute the Euler rotations; To install pyquaternion:

pip install pyquaternion

The program generate_fractal.py executes coord_generator.out and generates a file coords.data which contains the coordinates of the fractal monomers. Then, these coordinates are rotated by the Euler angles alpha, beta and gamma. Finally, the function gen_frac generates executes mstm.out using a file which includes all the angles that you want to calculate. By default the file is "scat_angles.csv". A file is created "mueller_log" in the save location (which by default is "test"). The file mueller_log and log contains the scattering phase function. 

To use the python function:
- mueller_matrix,coordinates = gen_frac(size_parameter,alpha,beta,gamma) 


IF you find this code useful, please cite:
- Piedra, Patricio, Aimable Kalume, Evgenij Zubko, Daniel Mackowski, Yong-Le Pan, and Gorden Videen. "Particle-shape classification using light scattering: An exercise in deep learning." Journal of Quantitative Spectroscopy and Radiative Transfer 231 (2019): 140-156.


