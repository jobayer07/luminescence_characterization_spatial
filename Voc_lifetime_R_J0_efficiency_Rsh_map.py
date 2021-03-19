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
import cv2

#Global values
V_01sun=0.494036
V_02sun=0.540099
V_1sun=0.607716
V_1sun_short=0

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

plt.figure(1)
im1 = plt.imshow(V_xy,vmin = 0.58, vmax = 0.63, cmap = 'inferno')
cbar = plt.colorbar()
plt.title('Voc Map')
plt.xlabel('Pixel Index-Horizontal')
plt.ylabel('Pixel Index-Vertical')


#------------------------Carrier lifetime Calculation------------------------------------

NA=1.49e16*1e6              #doping concentration (p-doped), /cm^3 -> /m^3
ni=1.45e10*1e6              #intrinsic carrier concentration/cm^3
Voc=V_xy
Jph=350                      #photocurrent=Jsc mA/cm^2 -> A/m^2
w=200e-6              #width of the cell, 200 um
q=1.6e-19             #electron charge
Vt=25.85e-3         #25mV=KT/q

Del_n=(np.sqrt(NA**2+4*ni**2*np.exp(Voc/Vt))-NA)/2
numerator=ni**2*np.exp(Voc/Vt)
denominator=Jph*(NA+Del_n)/(q*w)

Tau_xy=numerator/denominator

plt.figure(2)
im2 = plt.imshow(Tau_xy,vmin = 0.10e-6, vmax = 0.40e-4, cmap = 'inferno')
cbar = plt.colorbar()
plt.title('Lifetime Map')
plt.xlabel('Pixel Index-Horizontal')
plt.ylabel('Pixel Index-Vertical')



#--------------------------Series Resistance and Dark Saturation Current Density Calculation--------
Vapp1=V_1sun
Vapp2=V_1sun_short
V1_xy=V_xy

#V2_xy calculation
IL2=1
PL_image_2 = pd.read_csv('im_1_short.txt', header=None,delimiter = '\t')
PL_image_2[PL_image_2<0]=1
I_C_xy = 10*PL_image_2.as_matrix()                       #0.2sun
V2_xy=VT*np.log(abs((I_C_xy-B_xy*IL2)/C_xy))

#J0 calculation
factor=(Vapp1-V1_xy)/(Vapp2-V2_xy)
J0_xy=Jph*(factor-1)/(factor*np.exp(V2_xy/VT)-np.exp(V1_xy/VT))
J0_xy[J0_xy=='nan']=0

plt.figure(3)
im2 = plt.imshow(J0_xy, vmin = 1e-9, vmax = 5e-8, cmap = 'inferno')
cbar = plt.colorbar()
plt.title('Dark Saturation Current Density Map')
plt.xlabel('Pixel Index-Horizontal')
plt.ylabel('Pixel Index-Vertical')


#R_xy calculation
R_xy=(Vapp1-V1_xy)/(J0_xy*np.exp(V1_xy/VT)-Jph)
R_xy[R_xy>1e20]=0
R_xy[R_xy<0]=0

plt.figure(4)
im3 = plt.imshow(R_xy, vmin = 0.0013, vmax = 0.002, cmap = 'inferno')
cbar = plt.colorbar()
plt.title('Series Resistance Map')
plt.xlabel('Pixel Index-Horizontal')
plt.ylabel('Pixel Index-Vertical')

#--------------------------------------Efficiency Calculation--------------------------
NFF = pd.read_csv('NFF.txt', header=None,delimiter = '\t')
FF_xy = 10*PL_image_2.as_matrix()

Eff_xy=Jph*V_xy*100/1000
plt.figure(5)
im4 = plt.imshow(Eff_xy, vmin = 20, vmax = 22, cmap = 'inferno')
cbar = plt.colorbar()
plt.title('Efficiency Map')
plt.xlabel('Pixel Index-Horizontal')
plt.ylabel('Pixel Index-Vertical')

#--------------------------------------Efficiency Calculation--------------------------
del2_V2=cv2.Laplacian(V1_xy,cv2.CV_64F)
del2_V1=cv2.Laplacian(V2_xy,cv2.CV_64F)

fact2=del2_V1/del2_V2
J0_sh_xy=Jph*(fact2-1)/(fact2*np.exp(V2_xy/VT)-np.exp(V1_xy/VT))
J0_sh_xy[J0_sh_xy<-1]=0
Ro_sh_xy=del2_V1/(J0_sh_xy*np.exp(V1_xy/VT)-Jph)
Ro_sh_xy[Ro_sh_xy<-1]=0

plt.figure(6)
im5 = plt.imshow(J0_sh_xy, vmin = -0.05e-7, vmax = 5e-8, cmap = 'inferno')
cbar = plt.colorbar()
plt.title('Dark Saturation Current Density (model 2)')

plt.figure(7)
im6 = plt.imshow(Ro_sh_xy, vmin = -0.0001, vmax = 0.0001, cmap = 'inferno')
cbar = plt.colorbar()
plt.title('Shunt Resistance')