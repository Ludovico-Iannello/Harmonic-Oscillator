import pylab as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.special as ssp
plt.rc('font',size=16)

beta_omega=10.0

x1=np.loadtxt('/home/ludovico/Harmonic-Oscillator/data_varing_eta/fort.11', unpack=False)
x2=np.loadtxt('/home/ludovico/Harmonic-Oscillator/data_varing_eta/fort.66', unpack=False)
x3=np.loadtxt('/home/ludovico/Harmonic-Oscillator/data_varing_eta/fort.220', unpack=False)

a=[beta_omega/(len(x1)-1),beta_omega/(len(x2)-1),beta_omega/(len(x3)-1)]

n1=np.arange(0,len(x1))
n2=np.arange(0,len(x2))         
n3=np.arange(0,len(x3))

#path
plt.figure(1,figsize=[10,8])
plt.minorticks_on()
plt.title('path-integral')
plt.errorbar(n1*a[0],x1,label='N=10')
plt.errorbar(n2*a[1],x2,label='N=60')
plt.errorbar(n3*a[2],x3,label='N=200')
plt.xlabel(r'$\omega\tau$ ')
plt.ylabel('y ')
plt.legend()
plt.tight_layout()
'''
## waves function
path= np.loadtxt(r'/home/ludovico/Harmonic-Oscillator/meas_out')

N=200
#trasformo la matrice in un vettore dato che nono importa l'ordine
y = np.reshape(path, N*10**6)
def G(x, m):
    #Harmonic oscillator solution
    psi = (1/(np.pi)**(1/4))*(1/np.sqrt((2**m)*ssp.gamma(m+1)))*ssp.eval_hermite(m, x)*np.exp(-(x**2)/2)
    return psi
    
plt.figure(2,figsize=[10,8])
x = np.linspace(-3,3, 10000)
plt.hist(y, bins=1000, density=True)#, label=r'$\beta \omega$ = 20; $\eta$=0.1')
plt.plot(x, abs(G(x, 0))**2,lw=2.5)
plt.xlabel('y')
plt.ylabel(r'$|\psi_0|^2$')
plt.tight_layout()
'''




#plot varing eta
beta_omega=20.0
def fit (x,a,b):
    return a+b*x**2#+c*x**3
    
y2,var_y2,Dy2,var_Dy2=np.loadtxt('/home/ludovico/Harmonic-Oscillator/data_varing_eta_8/observables_all_paths', unpack=True)

eta=[]
N=[]
for i in range(1,len(y2)+1,1):
    N.append(i*beta_omega)
    eta.append(beta_omega/N[i-1])

eta=np.array(eta)
plt.figure(3,figsize=[10,8])

pars, covm = curve_fit(fit, eta, y2, sigma=var_y2, absolute_sigma=True)
err = np.sqrt(covm.diagonal())
chisq = np.sum(((y2 - fit(eta, *pars))/var_y2)**2)
ndof = len(eta) - len(pars)

xx=np.linspace(min(eta),max(eta),1000)
plt.minorticks_on()

plt.subplot(3,1,(1,2))
plt.title(r'fit $\langle y^2 \rangle$ al variare di $\eta$')
plt.ylabel(r'$y^2$ ')
plt.errorbar(eta,y2,var_y2,marker='.', linestyle='',capsize=4)
plt.plot(xx,fit(xx,*pars))
plt.subplot(3,1,3)
plt.axhline(0)
plt.errorbar(eta,((y2 - fit(eta, *pars))/var_y2),1,marker='.',capsize=4,linestyle='')
plt.xlabel(r'$\eta$')
plt.ylabel(r'res.norm')
plt.tight_layout()

print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
print('parametri fit: y^2 medio, A')
print(pars)
print('errori fit')
print(err)

#plot correction to K
plt.figure(4,figsize=[12,6])
plt.subplot(121)
plt.title(r'Termine cinetico non norm.')
plt.ylabel(r'$-\frac{\langle \Delta y^2 \rangle}{2\eta^2}$')
plt.xlabel(r'$\eta$')
plt.errorbar(eta, -Dy2/(2*eta**2), var_Dy2/(2*eta**2), fmt='.',capsize=4)

plt.subplot(122)
plt.title(r'Termine cinetico norm. e potenziale')
plt.xlabel(r'$\eta$')
plt.ylabel(r'$\frac{1}{2 \eta}-\frac{\langle \Delta y^2 \rangle}{2\eta^2}$')
plt.errorbar(eta, 1/(2*eta) - Dy2/(2*eta**2), var_Dy2/(2*eta**2), fmt='.',capsize=4,label='termine cinetico')
plt.errorbar(eta, y2/2,var_y2/2, fmt='.',capsize=4,label='termine potenziale')
plt.legend()
plt.tight_layout()

# plot varing beta
y2,var_y2,Dy2,var_Dy2=np.loadtxt('/home/ludovico/Harmonic-Oscillator/data_varing_beta/observables_all_paths_varing_beta', unpack=True)

eta=0.1
beta_omega=[]
N=[]
for i in range(1,21,1):
    beta_omega.append(i)
    N.append(beta_omega[i-1]/eta)


plt.figure(5,figsize=[10,8])

def fit (x,a):
    return a+1/(np.exp(x)-1)

U = 1/(2*eta) - Dy2/(2*eta**2) +y2/2
dU= np.sqrt((var_y2/2)**2 + (var_Dy2/(2*eta**2))**2)

pars, covm = curve_fit(fit, beta_omega, U, sigma=var_y2, absolute_sigma=True)
err = np.sqrt(covm.diagonal())
chisq = np.sum(((U - fit(beta_omega, *pars))/dU)**2)
ndof = len(beta_omega) - len(pars)

xx=np.linspace(min(beta_omega),max(beta_omega),1000)
plt.minorticks_on()

plt.subplot(3,1,(1,2))
plt.title(r'fit U al variare della temperatura')
plt.ylabel(r'U')
plt.errorbar(beta_omega,U,dU,marker='.', linestyle='',capsize=4)
plt.plot(xx,fit(xx,*pars))
plt.subplot(3,1,3)
plt.axhline(0)
plt.errorbar(beta_omega,((U - fit(beta_omega, *pars))/dU),1,marker='.',capsize=4,linestyle='')
plt.xlabel(r'$\beta \omega$')
plt.ylabel(r'res.norm')
plt.tight_layout()

print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
print('parametri fit: y^2 medio, A')
print(pars)
print('errori fit')
print(err)


uno_su_beta_omega=[]
for i in range(1,21,1):
    uno_su_beta_omega.append(1.0/i)

plt.figure(6,figsize=[10,8])

def fit (x,a):
    return a+1/(np.exp(x)-1)

U = 1/(2*eta) - Dy2/(2*eta**2) +y2/2
dU= np.sqrt((var_y2/2)**2 + (var_Dy2/(2*eta**2))**2)

U_diff=U-U[19]*np.ones(len(U))
dU_diff=np.sqrt(dU**2+dU[0]**2*np.ones(len(U)))


xx=np.linspace(min(uno_su_beta_omega),max(uno_su_beta_omega),1000)
plt.minorticks_on()


plt.ylabel(r'$U[T,\eta]-U[0,\eta]$')
plt.errorbar(uno_su_beta_omega,U_diff,dU_diff,marker='.', linestyle='',capsize=4)
plt.errorbar(xx,fit(1/xx,0),capsize=4)

plt.xlabel(r'$\frac{    1}{\beta \omega}$')
plt.tight_layout()


plt.show()
