import pylab as plt
import numpy as np
from scipy.optimize import curve_fit

y2,var_y2,Dy2,var_Dy2=np.loadtxt('/home/ludovico/Harmonic-Oscillator/data_varing_eta_2/observables_all_paths', unpack=True)

plt.rc('font',size=16)
misura=2

#plot varing eta
plt.figure(1,figsize=[10,8])


if misura==1:
    #data=[1,2,4,6,10,12,19]
    #data=[1,2,3,5,6,10,16]
    data=[1,2,3,4,5,6,12,16]
    beta_omega=20.0
else:
    data=[1,2,3,4,6,10,15]
    data=[1,2,3,4,5,6,10,15]
    beta_omega=10.0
eta=[]
N=[]
corr_matr=[]
dcorr_matr=[]
k_matr=[]
pars_arr=[]
err_arr=[]
eta2=[]

for i in data:
    print(i*10)
    N.append(i*beta_omega)
    eta.append(beta_omega/(i*beta_omega))
    if misura==1:
        k,corr,dcorr=np.loadtxt('/home/ludovico/Harmonic-Oscillator/data_varing_eta_7/fort.'+str(i*10), unpack=True)
    else:
        k,corr,dcorr=np.loadtxt('/home/ludovico/Harmonic-Oscillator/data_varing_eta_2/fort.'+str(i*10), unpack=True)
        corr=corr-(y2[i-1]**2)*np.ones(len(corr))
        #dcorr=np.sqrt(dcorr**2+var_y2[i-1]**2*np.ones(len(corr)))

    corr_matr.append(corr)
    dcorr_matr.append(dcorr)
    k_matr.append(k)

def fit (x,a,ene):
    return a*np.exp(-ene*x)

marker_styles = ['s', 'D', '^', 'v', '<', '>', 'p', '*', 'h', '+', 'x', '.', ',', '1', '2', '3', '4', '|']

for i in range(0,len(corr_matr),1):
    pars, covm = curve_fit(fit, k_matr[i]*beta_omega/N[i], corr_matr[i], sigma=dcorr_matr[i], absolute_sigma=True)
    err = np.sqrt(covm.diagonal())
    chisq = np.sum(((corr_matr[i] - fit(k_matr[i]*beta_omega/N[i], *pars))/dcorr_matr[i])**2)
    ndof = len(eta) - len(pars)
    if misura ==1:    
        xx=np.linspace(0,10,1000)
    else: 
        xx=np.linspace(0,5,1000)   
    plt.minorticks_on()
    
    
    plt.xlabel(r'$\omega\tau$')
    plt.errorbar(k_matr[i]*beta_omega/N[i],corr_matr[i],dcorr_matr[i],capsize=4,
                 marker=marker_styles[i],linestyle='',label=r'$\eta=$'+str((round(eta[i],2))))
    if misura==1:
        plt.errorbar(xx,0.5*np.exp(-xx),linestyle='--',color='black')
        plt.title(r'fit esponenziale di $\langle y(\tau)y(0) \rangle$')
        plt.ylabel(r'$\langle y(\tau)y(0) \rangle$')
    else: 
        plt.errorbar(xx,0.5*np.exp(-2*xx),linestyle='--',color='black')
        plt.title(r'fit esponenziale di $\langle y^2(\tau)y^2(0) \rangle$')
        plt.ylabel(r'$\langle y^2(\tau)y^2(0) \rangle$')
    plt.tight_layout()

    print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
    print('parametri fit: ')
    print(pars)
    print('errori fit')
    print(err)
    pars_arr.append(pars[1])
    err_arr.append(err[1])
    eta2.append(eta[i]**2)
    plt.legend()

eta2=np.array(eta2)
pars_arr=np.array(pars_arr)
err_arr=np.array(err_arr)
plt.figure(2,figsize=[10,8])
def fit (x,a,b):
    return a+b*x
pars, covm = curve_fit(fit, eta2, pars_arr, sigma=err_arr, absolute_sigma=True)
err = np.sqrt(covm.diagonal())
chisq = chisq = np.sum(((pars_arr - fit(eta2, *pars))/err_arr)**2)
ndof = len(eta2) - len(pars)
xx=np.linspace(min(eta2),max(eta2),1000)
plt.minorticks_on()

plt.subplot(3,1,(1,2))
if misura==1:
    plt.ylabel(r'$\frac{E_1-E_0}{\hbar\omega}$ ')
else:
    plt.ylabel(r'$\frac{E_2-E_0}{\hbar\omega}$ ')

plt.errorbar(eta2,pars_arr,err_arr,marker='.', linestyle='',capsize=4)
plt.plot(xx,fit(xx,*pars))
plt.subplot(3,1,3)
plt.axhline(0)
plt.errorbar(eta2,((pars_arr - fit(eta2, *pars))/err_arr),1,marker='.',capsize=4,linestyle='')
plt.xlabel(r'$\eta^2$')
plt.ylabel(r'res.norm')
plt.tight_layout()

print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
print('parametri fit: y^2 medio, A')
print(pars)
print('errori fit')
print(err)

plt.show()
