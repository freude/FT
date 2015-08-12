#!/usr/bin/python

import numpy as np

#-----------------------------constants-----------------------------------

h=1.05e-34;
q=1.6e-19;
m0=9.11e-31;
eps0=8.85e-12;

#----------------------silicon material parameters------------------------

mper=0.19*m0
mpar=0.98*m0
eps1=11.4
a_Si=0.541*1e-9

#--------------------------scaled atomic units----------------------------

k0 = 0.85*2*np.pi/a_Si
ab = 4*np.pi*h/q*h/q*eps1*eps0/mper
E_Har = q*q/ab/4/np.pi/eps0/eps1

