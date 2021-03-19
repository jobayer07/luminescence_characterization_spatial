"""@author: mhossain42
this is totally my approach
"""
#                           Code With Background 
#--------------------------Imports--------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#--------------------------Calibration-------------------------------------------
PL_image_1 = pd.read_csv('im_02.txt', header=None,delimiter = '\t')
I_B_xy = 2*PL_image_1.as_matrix()                         #0.1sun
PL_image_2 = pd.read_csv('im_01.txt', header=None,delimiter = '\t')
I_C_xy = PL_image_2.as_matrix()                       #0.2sun

Vb=0.540099
I_LB=0.1                                                    #0.1sun
Vc=0.494036                                                   #It varies with injection level
I_LC=0.2                                                    #0.2sun
VT=25.85e-03                                                #KT/q = 25meV

B_xy=(I_B_xy*np.exp(-Vb)-I_C_xy*np.exp(-Vc))/(I_LB*np.exp(-Vb)-I_LC*np.exp(-Vc))                                           #If we want to utilize short circit PL, Background not zero
#B_xy=0                                                      #Background=0
C_xy=(I_C_xy-B_xy*I_LC)/np.exp(Vc/VT)

#--------------------------Main Code------------------------------------------
current_data = pd.read_csv('im_1.txt', header=None,delimiter = '\t')   #sample: 1-100-011
I_xy = 10*current_data.as_matrix()
IL=1                                                       #text file 31517, 0.2sun

V_xy=VT*np.log((I_xy-B_xy*IL)/C_xy)
m=np.min(V_xy)
n=np.max(V_xy)

#im1=plt.imshow(I_xy, cmap = 'inferno')
im = plt.imshow(V_xy,vmin = 0.55, vmax = 0.58, cmap = 'inferno')
cbar = plt.colorbar()

