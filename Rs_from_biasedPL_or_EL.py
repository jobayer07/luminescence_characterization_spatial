"""
@author: Mohammad Jobayer Hossain
"""

#--------------------------Imports--------------------------------------------
#R=abs(V_d,xy-V_mpp)/J_mpp
#-----------------------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage

#------------------units-------------

fA=1e-15
nA=1e-9
mA=1e-3
VT=VT=25.85e-03

#-------------------------------------------

Voc_low_sun=0.554   #volt
Jsc=39.087609*mA    #A/cm2
Vmpp=517*mA
J01=0



mf=0.0942               #metal fraction
J01_pass_f=200*fA       #J01 passivated region, front  
J01_metal_f=6900*fA

J01_pass_r=0
J01_metal_r=500*fA

J01_pass=J01_pass_f+J01_pass_r
J01_metal=J01_metal_f+J01_metal_r

J01=J01_pass*(1-mf)+ J01_metal*mf

#----------------------------------Functions-----------------------------------------
def generate_Voltage_image(I_xy, IC_xy, Voc_low_sun):                                              #KT/q = 25meVn
    C_xy=(IC_xy)/np.exp(Voc_low_sun/VT)
    V_xy=VT*np.log(abs((I_xy)/C_xy))
    V_xy[V_xy > 100] = 0 #replacing inf with Voc, if v_xy>100
    V_xy[V_xy<0]=0
    return V_xy

#--------------------------Voc map first: Without Background--------------------------------------------
pl_folder='D:/STUDY_CREOL/Research with Davis/Characterization_Other_Projects/DOE_PVRD2_ContactDegradation/Rs_from_PL_EL/Part 1/Foshan_biased_PL'
el_folder='D:/STUDY_CREOL/Research with Davis/Characterization_Other_Projects/DOE_PVRD2_ContactDegradation/Rs_from_PL_EL/Part 1/Foshan_EL'

pl_1sun_open = pd.read_csv(str(pl_folder)+r'/840_Illum_1_Voc_617_txt.txt', header=None, delim_whitespace=True).to_numpy()
pl_0p1sun_open = pd.read_csv(str(pl_folder)+r'/840_Illum_0.1_Voc.txt', header=None, delim_whitespace=True).to_numpy()
pl_vmpp = pd.read_csv(str(pl_folder)+r'/840_Illum_1_Vmpp_524_txt.txt', header=None, delim_whitespace=True).to_numpy()
el_vmpp = pd.read_csv(str(el_folder)+r'/840_Illum_0_Vmpp_524_txt.txt', header=None, delim_whitespace=True).to_numpy()

def calculate_rs(luminescence_vmpp, sig_type):
    VT=25.85e-03 
    Vd_xy=generate_Voltage_image(luminescence_vmpp, pl_0p1sun_open, Voc_low_sun)      #Diode voltage
    
    if (sig_type=='pl'):
        Jmpp=Jsc-J01*(np.exp(Vmpp/VT)-1)
    elif (sig_type=='el'):
        Jmpp=J01*(np.exp(Vmpp/VT)-1)
    
    Rs=abs(Vd_xy-Vmpp)/Jmpp
    return Rs

R_xy=calculate_rs(el_vmpp, 'el')

plt.figure(1)
#plt.imshow(R_xy, vmin=0, vmax=3, cmap = 'inferno')
plt.imshow(R_xy, vmin=0, vmax=70, cmap = 'inferno')
cbar = plt.colorbar()
plt.title('Rs ($\Omega.cm^2$)')
plt.axis('off')