import numpy as np

h=6.62607015e-34 # J/Hz
c=299792458 # m/s
kB=1.380649e-23 # m^2 kg/s^2/K

def planck(lamb, T):
        return ((2*h*c**2)/(lamb**5))/(np.exp((h*c)/(lamb*kB*T))-1)

def planck_dT(lamb, T):
        return ((2 * h**2 * c**3)/(lamb**6 * kB * T**2)) * np.exp((h*c)/(lamb*kB*T)) * (np.exp((h*c)/(lamb*kB*T))-1)**(-2)

def bandpass_correction(lamb_lims, T):
        lamb=np.linspace(lamb_lims[0], lamb_lims[1], 1000)
        B=planck(lamb, T)

        # calculate \int B_\lambda d\lambda
        dlamb=lamb[1:]-lamb[:-1]
        med_B=(B[1:]+B[:-1])/2
        int_B=np.sum(dlamb*med_B)

        # calculate \int (\partial B_\lambda/ \partial T) d\lambda
        dBdT=planck_dT(lamb, T)
        med_dBdT=(dBdT[1:]+dBdT[:-1])/2
        int_dBdT=np.sum(dlamb*med_dBdT)

        bpc = (T/4)*(int_dBdT/int_B)

        return bpc
