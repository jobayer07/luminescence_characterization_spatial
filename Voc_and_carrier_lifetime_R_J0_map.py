"""
Open Circuit parameters
V(0.1 sun)=0.494036
V(0.2 sun)=0.540099
V(1 sun)=0.607716
"""

#--------------------------Imports--------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#Global values
V_01sun=0.494036
V_02sun=0.540099
V_1sun=0.607716

#--------------------------Voc map first: Without Background--------------------------------------------
PL_image_1 = pd.read_csv('im_01.txt', header=None,delimiter = '\t')
PL_image_1[PL_image_1<0]=1
I_C_xy = 1*PL_image_1.as_matrix()                       #0.1sun

Vc=0.494036                                                   #It varies with injection level
I_LC=0.1                                                    #0.2sun open
VT=25.85e-03                                                #KT/q = 25meV
B_xy=0                                                      #Background=0
C_xy=(I_C_xy-B_xy*I_LC)/np.exp(Vc/VT)

current_data = pd.read_csv('im_1.txt', header=None,delimiter = '\t')   #sample: 1-100-011
current_data[current_data<0]=1
I_xy = 10*current_data.as_matrix()
IL=1                                                       #text file 31517, 0.2sun


V_xy=VT*np.log(abs((I_xy-B_xy*IL)/C_xy))
V_xy[V_xy > 100] = 0.5                      #replacing inf with 0.5, v_xy=v_xy, if v_xy>100

#im = plt.imshow(V_xy,vmin = 0.53, vmax = 0.58, cmap = 'inferno')
#cbar = plt.colorbar()

#------------------------Carrier lifetime Calculation------------------------------------

NA=1.49e16*1e6              #doping concentration (p-doped), /cm^3 -> /m^3
ni=1.45e10*1e6              #intrinsic carrier concentration/cm^3
Voc=V_xy
Jph=35e-1                      #photocurrent=Jsc mA/cm^2 -> A/m^2
w=200e-6              #width of the cell, 200 um
q=1.6e-19             #electron charge
Vt=25.85e-3         #25mV=KT/q

Del_n=(np.sqrt(NA**2+4*ni**2*np.exp(Voc/Vt))-NA)/2
numerator=ni**2*np.exp(Voc/Vt)
denominator=Jph*(NA+Del_n)/(q*w)

Tau_xy=numerator/denominator
#im2 = plt.imshow(Tau_xy,vmin = 0.10e-2, vmax = 0.40e-2, cmap = 'inferno')
#cbar = plt.colorbar()

#--------------------------Series Resistance and Dark Saturation Current Density Calculation--------
Vapp1=V_1sun
Vapp2=V_02sun
V1_xy=V_xy

#V2_xy calculation
IL2=0.2
PL_image_2 = pd.read_csv('im_02.txt', header=None,delimiter = '\t')
PL_image_2[PL_image_2<0]=1
I_C_xy = 2*PL_image_2.as_matrix()                       #0.2sun
V2_xy=VT*np.log(abs((I_xy-B_xy*IL2)/C_xy))

#J0 calculation
factor=(Vapp1-V1_xy)/(Vapp2-V2_xy)
J0_xy=Jph*(factor-1)/(factor*np.exp(V2_xy/VT)-np.exp(V1_xy/VT))
J0_xy[J0_xy=='nan']=0
#im2 = plt.imshow(J0_xy, vmin = 1e-13, vmax = 5e-10, cmap = 'inferno')
#cbar = plt.colorbar()

#R_xy calculation
R_xy=(Vapp1-V1_xy)/(J0_xy*np.exp(V1_xy/VT)-Jph)
R_xy[R_xy>1e20]=0
R_xy[R_xy<0]=0
#im2 = plt.imshow(R_xy, vmin = 1e11, vmax = 7e13, cmap = 'inferno')
#cbar = plt.colorbar()

#--------------------------Series Resistance and Dark Saturation Current Density Calculation--------
'''

import numpy as np

#first of all calculate laplacian function del^2
def del2(M):
    dx = 1
    dy = 1
    rows, cols = M.shape
    dx = dx * np.ones ((1, cols - 1))
    dy = dy * np.ones ((rows-1, 1))

    mr, mc = M.shape
    D = np.zeros ((mr, mc))

    if (mr >= 3):
        ## x direction
        ## left and right boundary
        D[:, 0] = (M[:, 0] - 2 * M[:, 1] + M[:, 2]) / (dx[:,0] * dx[:,1])
        D[:, mc-1] = (M[:, mc - 3] - 2 * M[:, mc - 2] + M[:, mc-1]) \
            / (dx[:,mc - 3] * dx[:,mc - 2])

        ## interior points
        tmp1 = D[:, 1:mc - 1] 
        tmp2 = (M[:, 2:mc] - 2 * M[:, 1:mc - 1] + M[:, 0:mc - 2])
        tmp3 = np.kron (dx[:,0:mc -2] * dx[:,1:mc - 1], np.ones ((mr, 1)))
        D[:, 1:mc - 1] = tmp1 + tmp2 / tmp3

    if (mr >= 3):
        ## y direction
        ## top and bottom boundary
        D[0, :] = D[0,:]  + \
            (M[0, :] - 2 * M[1, :] + M[2, :] ) / (dy[0,:] * dy[1,:])

        D[mr-1, :] = D[mr-1, :] \
            + (M[mr-3,:] - 2 * M[mr-2, :] + M[mr-1, :]) \
            / (dy[mr-3,:] * dx[:,mr-2])

        ## interior points
        tmp1 = D[1:mr-1, :] 
        tmp2 = (M[2:mr, :] - 2 * M[1:mr - 1, :] + M[0:mr-2, :])
        tmp3 = np.kron (dy[0:mr-2,:] * dy[1:mr-1,:], np.ones ((1, mc)))
        D[1:mr-1, :] = tmp1 + tmp2 / tmp3

    return D / 4

#Main calculations for J0 starts
factor2=del2(V1_xy)/del2(V2_xy)
J0_2_xy=Jph*(factor2-1)/(factor2*np.exp(V2_xy/VT)-np.exp(V1_xy/VT))

im2 = plt.imshow(factor2, vmin = 0, vmax = 1, cmap = 'inferno')
cbar = plt.colorbar()
'''

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