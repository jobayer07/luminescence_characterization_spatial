#--------------------------Imports--------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

#import torch

Voc_low_sun=0.588   #volt
Jsc=38.020701 #mA/cm2
Jmpp=35.639156

#pl0_folder='D:/STUDY_CREOL/Research with Davis/Characterization_Other_Projects/DOE_PVRD2_ContactDegradation/Rs_from_PL_EL/Part 1/Foshan_biased_PL/'

pl_folder='D:/STUDY_CREOL/Research with Davis/Characterization_Other_Projects/DOE_PVRD2_ContactDegradation/Rs_from_PL_EL/Cell_with_defects/M9004C002101/Biased PL/'
el_folder='D:/STUDY_CREOL/Research with Davis/Characterization_Other_Projects/DOE_PVRD2_ContactDegradation/Rs_from_PL_EL/Cell_with_defects/M9004C002101/EL/'
pl_0p1sun_open = pd.read_csv(str(pl_folder)+r'M9004C002101_0.1.txt', header=None, delim_whitespace=True).to_numpy()
#----------------------------------Functions-----------------------------------------
def generate_Voltage_image(I_xy):
    VT=25.85e-03                                                #KT/q = 25meVn
    
    IC_xy=pl_0p1sun_open
    C_xy=(IC_xy)/np.exp(Voc_low_sun/VT)
    
    V_xy=VT*np.log(abs((I_xy)/C_xy))
    V_xy[V_xy > 100] = 0 #replacing inf with Voc, if v_xy>100
    V_xy[V_xy<0]=0
    return V_xy

def interpolation(x,y, xnew): 
    tck = interpolate.splrep(x, y, s=0)  #s=0 -> No smoothing/regression required      
    ynew = interpolate.splev(xnew, tck, der=0)
    return ynew

def calculate_rs_image(i,j):
    pl_J=np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])*0.1*Jsc
    pl_V=np.array([plv1[i,j], plv2[i,j], plv3[i,j], plv4[i,j], plv5[i,j], plv6[i,j], plv7[i,j], plv8[i,j], plv9[i,j], plv10[i,j], plv11[i,j], plv12[i,j]])
    #---------------------------------------------------
    
    dark_JV=pd.read_table('D:/STUDY_CREOL/Research with Davis/Characterization_Other_Projects/DOE_PVRD2_ContactDegradation/Rs_from_PL_EL/Cell_with_defects/Dark IV/M9004C002101.txt', skiprows=2).to_numpy() #skipping first 1 row
    dark_V=dark_JV[:, 0]*1e-3
    dark_J=-dark_JV[:, 1]
    
    el_V=np.array([elv1[i,j], elv2[i,j], elv3[i,j], elv4[i,j], elv5[i,j], elv6[i,j], elv7[i,j], elv8[i,j], elv9[i,j], elv10[i,j], elv11[i,j], elv12[i,j], elv13[i,j], elv14[i,j], elv15[i,j], elv16[i,j], elv17[i,j], elv18[i,j], elv19[i,j], elv20[i,j], elv21[i,j]])
    V_applied=np.arange(0.50, 0.71, 0.01)
    el_J=interpolation(dark_V,dark_J, V_applied)   #this is actually J_applied
    
    #Rs calculation at Jmpp
    Vmpp_pl=interpolation(pl_J, pl_V, Jmpp)
    Vmpp_el=interpolation(el_J,el_V, Jmpp)
    Vmpp_darkIV=np.interp(Jmpp, dark_J, dark_V)
    Rs_pl_el_Jmpp=abs(Vmpp_pl-Vmpp_el)/(Jmpp*1e-3)    # Rs is in ohm-cm2 now
    Rs_el_darkIV_Jmpp=abs(Vmpp_darkIV-Vmpp_el)/(Jmpp*1e-3)
    
    #Rs calculation at Jsc
    V_pl_Jsc=interpolation(pl_J, pl_V, Jsc)
    V_el_Jsc=interpolation(el_J,el_V, Jsc)
    V_darkIV_Jsc=np.interp(Jsc, dark_J, dark_V)   # we calculate at Jsc
    Rs_pl_el_Jsc=abs(V_pl_Jsc-V_el_Jsc)/(Jsc*1e-3)        # Rs is in ohm-cm2 now
    Rs_el_darkIV_Jsc=abs(V_darkIV_Jsc-V_el_Jsc)/(Jsc*1e-3)

    #---------------------Plotting -----------------------------------
    '''
    plt.figure(1)
    plt.plot(el_V, el_J, '-r', label='EL converted JV')
    plt.plot(dark_V,dark_J, '-k', label='Dark JV (global)')
    plt.plot(pl_V, pl_J, label='PL converted JV')
    plt.legend()
    plt.xlabel('Voltage ($V$)', fontsize=12)
    plt.ylabel('Current Density ($mA/cm^2$)', fontsize=12)
    plt.ylim(0, 40)
    plt.xlim(0.45, 0.68)
    plt.title('At pixel ('+str(i)+','+str(j)+')')
    #Annotation
    plt.scatter(Vmpp_pl,Jmpp, cmap="hot", edgecolor='k', s = 20)
    plt.scatter(Vmpp_el,Jmpp, cmap="hot", edgecolor='k', s = 20)
    '''
    return Rs_pl_el_Jmpp, Rs_el_darkIV_Jmpp, Rs_pl_el_Jsc, Rs_el_darkIV_Jsc

def stat(var):
    mn=np.nansum(var)/np.count_nonzero(var)                 #only the nonzero values, excludes nans
    sd=np.nanstd(np.where(np.isclose(var,0), np.nan, var))  #only the nonzero values
    return mn,sd

#--------------------------Voc map first: Without Background--------------------------------------------

pl1 = pd.read_csv(str(pl_folder)+r'M9004C002101_0.1.txt', header=None, delim_whitespace=True).to_numpy()
pl2 = pd.read_csv(str(pl_folder)+r'M9004C002101_0.2.txt', header=None, delim_whitespace=True).to_numpy()
pl3 = pd.read_csv(str(pl_folder)+r'M9004C002101_0.3.txt', header=None, delim_whitespace=True).to_numpy()
pl4 = pd.read_csv(str(pl_folder)+r'M9004C002101_0.4.txt', header=None, delim_whitespace=True).to_numpy()
pl5 = pd.read_csv(str(pl_folder)+r'M9004C002101_0.5.txt', header=None, delim_whitespace=True).to_numpy()
pl6 = pd.read_csv(str(pl_folder)+r'M9004C002101_0.6.txt', header=None, delim_whitespace=True).to_numpy()
pl7 = pd.read_csv(str(pl_folder)+r'M9004C002101_0.7.txt', header=None, delim_whitespace=True).to_numpy()
pl8 = pd.read_csv(str(pl_folder)+r'M9004C002101_0.8.txt', header=None, delim_whitespace=True).to_numpy()
pl9 = pd.read_csv(str(pl_folder)+r'M9004C002101_0.9.txt', header=None, delim_whitespace=True).to_numpy()
pl10 = pd.read_csv(str(pl_folder)+r'M9004C002101_1.txt', header=None, delim_whitespace=True).to_numpy()
pl11 = pd.read_csv(str(pl_folder)+r'M9004C002101_1.1.txt', header=None, delim_whitespace=True).to_numpy()
pl12 = pd.read_csv(str(pl_folder)+r'M9004C002101_1.2.txt', header=None, delim_whitespace=True).to_numpy()

el1 = pd.read_csv(str(el_folder)+r'M9004C002101_0_500.txt', header=None, delim_whitespace=True).to_numpy()
el2 = pd.read_csv(str(el_folder)+r'M9004C002101_0_510.txt', header=None, delim_whitespace=True).to_numpy()
el3 = pd.read_csv(str(el_folder)+r'M9004C002101_0_520.txt', header=None, delim_whitespace=True).to_numpy()
el4 = pd.read_csv(str(el_folder)+r'M9004C002101_0_530.txt', header=None, delim_whitespace=True).to_numpy()
el5 = pd.read_csv(str(el_folder)+r'M9004C002101_0_540.txt', header=None, delim_whitespace=True).to_numpy()
el6 = pd.read_csv(str(el_folder)+r'M9004C002101_0_550.txt', header=None, delim_whitespace=True).to_numpy()
el7 = pd.read_csv(str(el_folder)+r'M9004C002101_0_560.txt', header=None, delim_whitespace=True).to_numpy()
el8 = pd.read_csv(str(el_folder)+r'M9004C002101_0_570.txt', header=None, delim_whitespace=True).to_numpy()
el9 = pd.read_csv(str(el_folder)+r'M9004C002101_0_580.txt', header=None, delim_whitespace=True).to_numpy()
el10 = pd.read_csv(str(el_folder)+r'M9004C002101_0_590.txt', header=None, delim_whitespace=True).to_numpy()
el11 = pd.read_csv(str(el_folder)+r'M9004C002101_0_600.txt', header=None, delim_whitespace=True).to_numpy()
el12 = pd.read_csv(str(el_folder)+r'M9004C002101_0_610.txt', header=None, delim_whitespace=True).to_numpy()
el13 = pd.read_csv(str(el_folder)+r'M9004C002101_0_620.txt', header=None, delim_whitespace=True).to_numpy()
el14 = pd.read_csv(str(el_folder)+r'M9004C002101_0_630.txt', header=None, delim_whitespace=True).to_numpy()
el15 = pd.read_csv(str(el_folder)+r'M9004C002101_0_640.txt', header=None, delim_whitespace=True).to_numpy()
el16 = pd.read_csv(str(el_folder)+r'M9004C002101_0_650.txt', header=None, delim_whitespace=True).to_numpy()
el17 = pd.read_csv(str(el_folder)+r'M9004C002101_0_660.txt', header=None, delim_whitespace=True).to_numpy()
el18 = pd.read_csv(str(el_folder)+r'M9004C002101_0_670.txt', header=None, delim_whitespace=True).to_numpy()
el19 = pd.read_csv(str(el_folder)+r'M9004C002101_0_680.txt', header=None, delim_whitespace=True).to_numpy()
el20 = pd.read_csv(str(el_folder)+r'M9004C002101_0_690.txt', header=None, delim_whitespace=True).to_numpy()
el21 = pd.read_csv(str(el_folder)+r'M9004C002101_0_700.txt', header=None, delim_whitespace=True).to_numpy()
#el22 = pd.read_csv(str(el_folder)+r'790_Illum_0_V_700_txt.txt', header=None, delim_whitespace=True).to_numpy()

plv1=generate_Voltage_image(pl1)
plv2=generate_Voltage_image(pl2)
plv3=generate_Voltage_image(pl3)
plv4=generate_Voltage_image(pl4)
plv5=generate_Voltage_image(pl5)
plv6=generate_Voltage_image(pl6)
plv7=generate_Voltage_image(pl7)
plv8=generate_Voltage_image(pl8)
plv9=generate_Voltage_image(pl9)
plv10=generate_Voltage_image(pl10)
plv11=generate_Voltage_image(pl11)
plv12=generate_Voltage_image(pl12)

elv1=generate_Voltage_image(el1)
elv2=generate_Voltage_image(el2)
elv3=generate_Voltage_image(el3)
elv4=generate_Voltage_image(el4)
elv5=generate_Voltage_image(el5)
elv6=generate_Voltage_image(el6)
elv7=generate_Voltage_image(el7)
elv8=generate_Voltage_image(el8)
elv9=generate_Voltage_image(el9)
elv10=generate_Voltage_image(el10)
elv11=generate_Voltage_image(el11)
elv12=generate_Voltage_image(el12)
elv13=generate_Voltage_image(el13)
elv14=generate_Voltage_image(el14)
elv15=generate_Voltage_image(el15)
elv16=generate_Voltage_image(el16)
elv17=generate_Voltage_image(el17)
elv18=generate_Voltage_image(el18)
elv19=generate_Voltage_image(el19)
elv20=generate_Voltage_image(el20)
elv21=generate_Voltage_image(el21)
#elv22=generate_Voltage_image(el22)

Rs_pl_el_Jmpp, Rs_el_darkIV_Jmpp, Rs_pl_el_Jsc, Rs_el_darkIV_Jsc=calculate_rs_image(200,215)

s,t=np.shape(el1)
Rs_pl_el_Jmpp=np.zeros([s,t])
Rs_el_darkIV_Jmpp=np.zeros([s,t])
Rs_pl_el_Jsc=np.zeros([s,t])
Rs_el_darkIV_Jsc=np.zeros([s,t])

for m in range (1,s-1):
    for n in range (1,t-1):
        Rs_pl_el_Jmpp[m,n], Rs_el_darkIV_Jmpp[m,n], Rs_pl_el_Jsc[m,n], Rs_el_darkIV_Jsc[m,n]=calculate_rs_image(m,n)



def mask_creation():
    mas=np.ones([s,t], dtype=float)
    i1=70
    i2=80
    i3=220
    i4=230
    i5=370
    i6=380
    i7=520
    i8=530
    for j in range(i1, i2):
        for i in range(s):
            mas[i,j]=0
    for j in range(i3, i4):
        for i in range(s):
            mas[i,j]=0
    for j in range(i5, i6):
        for i in range(s):
            mas[i,j]=0
    for j in range(i7, i8):
        for i in range(s):
            mas[i,j]=0
    return mas

mask=mask_creation()

#-----------------Image plotting-------------------------------
plt.figure(2)
plt.imshow(np.nan_to_num(Rs_pl_el_Jmpp)*mask, vmin=0, vmax=0.6, cmap='inferno')
plt.title('Rs_pl_el @ Jmpp (790°C)')
cbar = plt.colorbar()
plt.title('Rs_pl_el @ Jmpp ($\Omega.cm^2$)')
plt.axis('off') 

plt.figure(3)
plt.imshow(np.nan_to_num(Rs_el_darkIV_Jmpp)*mask, vmin=0, vmax=1.8, cmap='inferno')
cbar = plt.colorbar()
plt.title('Rs_el_darkIV @ Jmpp ($\Omega.cm^2$)')
plt.axis('off')     

plt.figure(4)
plt.imshow(np.nan_to_num(Rs_pl_el_Jsc)*mask, vmin=0, vmax=0.6, cmap='inferno')
cbar = plt.colorbar()
plt.title('Rs_pl_el @ Jsc ($\Omega.cm^2$)')
plt.axis('off')    
               

plt.figure(5)
plt.imshow(np.nan_to_num(Rs_el_darkIV_Jsc)*mask, vmin=0, vmax=1.8, cmap='inferno')
cbar = plt.colorbar()
plt.title('Rs_el_darkIV @ Jsc ($\Omega.cm^2$)')
plt.axis('off')

#--------------------------------Post Processing---------------------
mn_pl_el_Jmpp, sd_pl_el_Jmpp=stat(Rs_pl_el_Jmpp*mask)
mn_el_darkIV_Jmpp, sd_el_darkIV_Jmpp=stat(Rs_el_darkIV_Jmpp*mask)
mn_pl_el_Jsc, sd_pl_el_Jsc=stat(Rs_pl_el_Jsc*mask)
mn_el_darkIV_Jsc, sd_el_darkIV_Jsc=stat(Rs_el_darkIV_Jsc*mask)

all_stat=[[mn_pl_el_Jmpp, sd_pl_el_Jmpp], [mn_pl_el_Jsc, sd_pl_el_Jsc], [mn_el_darkIV_Jmpp, sd_el_darkIV_Jmpp], [mn_el_darkIV_Jsc, sd_el_darkIV_Jsc]]

