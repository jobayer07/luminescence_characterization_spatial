"""
Open Circuit parameters
V(0.1 sun)=0.494036
V(0.2 sun)=0.540099
V(1 sun)=0.607716
"""

#--------------------------Voc map first--------------------------------------------
#                           Code With Background 
#--------------------------Imports--------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#--------------------------Calibration-------------------------------------------
PL_image_1 = pd.read_csv('im_01_short.txt', header=None,delimiter = '\t')
PL_image_1[PL_image_1<0]=1
I_B_xy = PL_image_1.as_matrix()                         #0.1sun
PL_image_2 = pd.read_csv('im_01.txt', header=None,delimiter = '\t')
PL_image_2[PL_image_1<0]=1
I_C_xy = 1*PL_image_2.as_matrix()                       #0.2sun

I_LB=0.1                                                    #0.1sun short
Vc=0.494036                                                   #It varies with injection level
I_LC=0.1                                                    #0.2sun open
VT=25.85e-03                                                #KT/q = 25meV

B_xy=I_B_xy/I_LB                                           #If we want to utilize short circit PL, Background not zero
B_xy=0                                                      #Background=0
C_xy=(I_C_xy-B_xy*I_LC)/np.exp(Vc/VT)

#--------------------------Main Code------------------------------------------
current_data = pd.read_csv('im_1.txt', header=None,delimiter = '\t')   #sample: 1-100-011
current_data[current_data<0]=1
I_xy = 10*current_data.as_matrix()
IL=1                                                       #text file 31517, 0.2sun


V_xy=VT*np.log(abs((I_xy-B_xy*IL)/C_xy))
V_xy[V_xy > 100] = 0.5                      #replacing inf with 0.5, v_xy=v_xy, if v_xy>100
m=np.min(V_xy)
n=np.max(V_xy)

#im = plt.imshow(V_xy,vmin = 0.53, vmax = 0.58, cmap = 'inferno')
#cbar = plt.colorbar()

#------------------------Now Carrier lifetime-------------------------
NA=1.49e16*1e6              #doping concentration (p-doped), /cm^3 -> /m^3
ni=1.45e10*1e6              #intrinsic carrier concentration/cm^3
Voc=V_xy
Jph=35e-1                      #photocurrent=Jsc mA/cm^2 -> A/m^2
w=200e-6              #width of the cell, 200 um
q=1.6e-19             #electron charge
#--------------------------Main Code--------------------------------------------
Vt=25.85e-3         #25mV=KT/q
Del_n=(np.sqrt(NA**2+4*ni**2*np.exp(Voc/Vt))-NA)/2

numerator=ni**2*np.exp(Voc/Vt)
denominator=Jph*(NA+Del_n)/(q*w)

Tau_xy=numerator/denominator
im2 = plt.imshow(Tau_xy,vmin = 0.10e-2, vmax = 0.40e-2, cmap = 'inferno')
cbar = plt.colorbar()


'''
#Filter
from scipy.signal import butter, lfilter
fs = 1E9 # 1 ns -> 1 GHz
cutoff = 10E7 # 10 MHz
B, A = butter(1, cutoff / (fs / 2), btype='low') # 1st order Butterworth low-pass
filtered_signal = lfilter(B, A, R_xy, axis=0)

im2 = plt.imshow(filtered_signal, vmin = 0.16e-13, vmax = 0.16e-12, cmap = 'inferno')
cbar = plt.colorbar()
'''