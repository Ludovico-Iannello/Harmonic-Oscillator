import pylab as plt
import numpy as np
from scipy.optimize import curve_fit

plt.rc('font',size=16)


x1=np.loadtxt('fort.10', unpack=False)
x2=np.loadtxt('fort.20', unpack=False)
x3=np.loadtxt('fort.30', unpack=False)

n1=np.arange(1,len(x1)+1)
n2=np.arange(1,len(x2)+1)*10/(len(x2))
n3=np.arange(1,len(x3)+1)*10/(len(x3))
#path
plt.figure(1)
plt.minorticks_on()
plt.title('path-integral')
plt.errorbar(n1,x1)
plt.errorbar(n2,x2)
plt.errorbar(n3,x3)
plt.xlabel(r'$\tau$ [a.u.]')
plt.ylabel('y [a.u.]')
plt.show()
