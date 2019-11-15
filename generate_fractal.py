#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 13:56:24 2019

@author: drpiedra
"""

import pandas as pd
import numpy as np
from math import pi
from subprocess import STDOUT, check_output
import random as rm
from pyquaternion import Quaternion as qt


##### rotates a three dimensional array of vectos 
##### along angle alpha, beta and gamma
##### according to the definition of Euler angles
def rotate_euler(vector,alpha,beta,gamma):
    #rotates 3d vector in euler
    vector=np.array([vector.x.values,vector.y.values, vector.z.values]).T
    const=2*pi/360.

    xx=np.array([1,0,0])
    yy=np.array([0,1,0])
    zz=np.array([0,0,1])
    rotate_alpha=qt(axis=zz,angle=alpha*const)
    yy_2=rotate_alpha.rotate(yy)
    zz_2=rotate_alpha.rotate(zz)
    rotate_beta=qt(axis=yy_2,angle=beta*const)
    zz_3=rotate_beta.rotate(zz_2)
    rotate_gamma=qt(axis=zz_3,angle=gamma*const)

    rot_vec=[]
    for vec in vector:
        rotated=rotate_alpha.rotate(vec)
        rotated=rotate_beta.rotate(rotated)
        rotated=rotate_gamma.rotate(rotated)
        rot_vec.append(rotated)

    vector=np.array(rot_vec)
    return vector


#generates fractal scattering calculation for mstm
#it requires a file of theta,phi angles stored in scattering_angle_file 'scat_angles.csv'
#saves the mueller matrix in "save directory" =  test
def gen_frac(x,n_spheres,alpha,beta,gamma):
    #x,n_sphere,alpha,beta,gamma=5,100,0,0,0
    col_names=['theta','phi','s11','s12','s13','s14',\
                   's21','s22','s23','s24',\
                   's31','s32','s33','s34',\
                   's41','s42','s43','s44']

    #n_spheres=rm.randint(10,100) #number of spheres
    D=rm.uniform(1.7,2.2) #fractal dimension
    k0=(5/3.0)**(D/2.0) #constant

    save_directory='./test'
    #generate fractal geometry and save in file coords.dat, rotate and save in correct format
    cmd="printf %d,100,%.2f,%.2f,0,coords.dat | ./coord_generator.out"%(n_spheres,k0,D)
    output = check_output(cmd, stderr=STDOUT, timeout=10,shell=True)
    coords=pd.read_csv('coords.dat',names=['r','x','y','z'],delim_whitespace=True)
    rotated=rotate_euler(coords,alpha,beta,gamma)
    rotated=pd.DataFrame(rotated,columns=['x','y','z'])
    rotated['r']=coords.r
    rotated=rotated[['r','x','y','z']]
    rotated.to_csv('coords.dat',index=False,sep='\t',header=False)
    rotated[['x','y','z']].to_csv('%s/geom.csv'%save_directory,index=False) #save to read

    ##execute scattering
    f=open('scat_angles_input.inp','r')
    contents=f.read()
    contents= contents.split('\n')
    f.close()
    lscale=x*n_spheres**(-1./3.)
    contents[1]='%d'%n_spheres
    contents[7]='%.3fd0'%lscale
    out_file=save_directory+'/mueller_log'
    contents[9]=out_file
    
    #write a new input file
    g=open('input.inp','w')
    for line in contents:
        g.writelines(line+'\n')
    g.close()
    cmd2='./mstm.exe input.inp'
    output = check_output(cmd2, stderr=STDOUT, timeout=20,shell=True)
    mueller_matrix=pd.read_csv(out_file,compression=None,delim_whitespace=True, header=47)
    mueller_matrix.columns=col_names
    return mueller_matrix,coords[['x','y','z']]