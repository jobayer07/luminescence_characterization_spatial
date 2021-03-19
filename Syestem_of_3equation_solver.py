#--------------------------Imports--------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage  #used for image enlargement
from scipy import optimize

#---------------------------Solver--------------------------------------
m=np.matrix([[1,2,3],[4,5,6]])
x=m.flatten()
y=np.reshape(x, (2,3))
a1=1
b1=1
c1=3
a2=2
def myFunction(k):
        x=k[0]
        y=k[1]
        z=k[2]        
        F=np.empty((3))        
        F[0]=a1*x+b1*np.exp(y)+c1*z-17.38
        F[1]=a2*x+y+z+y*z-13
        F[2]=x**2+3*y+z**2-16
        return F
    

kGuess=np.array([15,20,15])
    
k=optimize.fsolve(myFunction,kGuess)
