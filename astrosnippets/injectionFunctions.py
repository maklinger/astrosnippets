import numpy as np

def powerlaw_array(Es, Emin, Emax, p, norm):
    powerlaw = (Es/Emin)**(-p) * np.exp(-(Es/Emax)) * np.greater_equal(Es, Emin)
    integral = np.log(Es[1]/Es[0]) * np.sum(Es**2 * powerlaw)
    return norm/integral * powerlaw

def thermal_array(Es, Eth, norm):
    powerlaw = (Es/Eth)**2 * np.exp(-(Es/Eth)) 
    integral = np.log(Es[1]/Es[0]) * np.sum(Es**2 * powerlaw)
    return norm/integral * powerlaw

def combined_thermal_powerlaw_array(Es, Eth, eps_th, Emax, p, eps_nth, norm, fac_PL_start=3):
    return thermal_array(Es, Eth, eps_th*norm) + powerlaw_array(Es, fac_PL_start*Eth, Emax, p, eps_nth*norm)