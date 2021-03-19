import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

temp='765'
Jsc=39.087600*1e-3 #A/cm2
voc_0p1sun=0.560
applied_voltage=np.array([0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60])


UT=25.85e-03
el_folder='D:/STUDY_CREOL/Research with Davis/Characterization_Other_Projects/DOE_PVRD2_ContactDegradation/Rs_from_PL_EL/Hinken_method/'+str(temp)+'_EL/'
pl_folder='D:/STUDY_CREOL/Research with Davis/Characterization_Other_Projects/DOE_PVRD2_ContactDegradation/Rs_from_PL_EL/Hinken_method/PL_folder/'

el1 = pd.read_csv(str(el_folder)+str(temp)+r'_Illum_0_V_490_txt.txt', header=None, delim_whitespace=True).to_numpy()
el2 = pd.read_csv(str(el_folder)+str(temp)+r'_Illum_0_V_500_txt.txt', header=None, delim_whitespace=True).to_numpy()
el3 = pd.read_csv(str(el_folder)+str(temp)+r'_Illum_0_V_510_txt.txt', header=None, delim_whitespace=True).to_numpy()
el4 = pd.read_csv(str(el_folder)+str(temp)+r'_Illum_0_V_520_txt.txt', header=None, delim_whitespace=True).to_numpy()
el5 = pd.read_csv(str(el_folder)+str(temp)+r'_Illum_0_V_530_txt.txt', header=None, delim_whitespace=True).to_numpy()
el6 = pd.read_csv(str(el_folder)+str(temp)+r'_Illum_0_V_540_txt.txt', header=None, delim_whitespace=True).to_numpy()
el7 = pd.read_csv(str(el_folder)+str(temp)+r'_Illum_0_V_550_txt.txt', header=None, delim_whitespace=True).to_numpy()
el8 = pd.read_csv(str(el_folder)+str(temp)+r'_Illum_0_V_560_txt.txt', header=None, delim_whitespace=True).to_numpy()
el9 = pd.read_csv(str(el_folder)+str(temp)+r'_Illum_0_V_570_txt.txt', header=None, delim_whitespace=True).to_numpy()
el10 = pd.read_csv(str(el_folder)+str(temp)+r'_Illum_0_V_580_txt.txt', header=None, delim_whitespace=True).to_numpy()
el11 = pd.read_csv(str(el_folder)+str(temp)+r'_Illum_0_V_590_txt.txt', header=None, delim_whitespace=True).to_numpy()
el12 = pd.read_csv(str(el_folder)+str(temp)+r'_Illum_0_V_600_txt.txt', header=None, delim_whitespace=True).to_numpy()

pl_1sun_open=pd.read_csv(str(pl_folder)+str(temp)+r'_Illum_1_Voc.txt', header=None, delim_whitespace=True).to_numpy()
pl_0p1sun_open=pd.read_csv(str(pl_folder)+str(temp)+r'_Illum_0.1_Voc.txt', header=None, delim_whitespace=True).to_numpy()

#----------------------------------------Function_one_pixel--------------------------
def generate_Voltage_image(I_xy, IC_xy, Voc_low_sun):
    VT=25.85e-03                                                #KT/q = 25meVn
    C_xy=(IC_xy)/np.exp(Voc_low_sun/VT)
    V_xy=VT*np.log(abs((I_xy)/C_xy))
    V_xy[V_xy > 100] = 0 #replacing inf with Voc, if v_xy>100
    V_xy[V_xy<0]=0
    return C_xy, V_xy

def calculate_mb(i,j):
    elx=np.array([el1[i,j], el2[i,j], el3[i,j], el4[i,j], el5[i,j], el6[i,j], el7[i,j], el8[i,j], el9[i,j], el10[i,j], el11[i,j], el12[i,j]])
    V=applied_voltage
    
    el_p=np.diff(elx)/np.diff(V)
    el=(elx[1:12]+elx[0:11])/2
    if (np.sum(el)==0):
        mb=100
    else:
        y=UT*el_p/el
        x=el_p        
        # Slope and intercept
        m,b=np.polyfit(x,y, deg=1)    #numpy.polyfit(x, y). m,b=solpe, intercept
        mb=m/b
        '''
        plt.figure(1)
        plt.plot(x, y, label='At ('+str(i)+',' +str(j)+ ')')
        plt.xlabel( r'$\phi^1$', fontsize=12)
        plt.ylabel( r'UT*$\phi^1/\phi$', fontsize=12)
        plt.legend()
        '''
    return mb

#a=calculate_mb(120,200)
#print(a)

#---------------------------------------Function_image--------------------------

I,J=np.shape(el1)
mb_xy=np.zeros([I,J])

for m in range (1,I-1):
    for n in range (1,J-1):
        mb_xy[m,n]=calculate_mb(m,n)

C_xy, Voc_xy=generate_Voltage_image(pl_1sun_open, pl_0p1sun_open, voc_0p1sun)
J0_xy=Jsc/(np.exp(Voc_xy/UT)-1)

Rs_xy=-(C_xy*mb_xy/J0_xy)
     

plt.figure(2)
plt.imshow(el12, cmap = 'inferno')
cbar = plt.colorbar()
plt.axis('off')

plt.figure(3)
plt.imshow(Rs_xy, vmin=0, vmax=2.5, cmap = 'inferno')
cbar = plt.colorbar()
plt.axis('off')

plt.figure(4)  
plt.hist(Rs_xy.ravel(), bins=100, range=(0,2.5), color='crimson')
#plt.title('Jsc histogram')

#calculate_rs_image(120,200)
#calculate_rs_image(300,300)
#calculate_rs_image(402,300)