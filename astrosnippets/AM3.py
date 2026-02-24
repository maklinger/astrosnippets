from .physicsConsts import eV2erg

def get_AM3_energy_grids(am3):
    Ee_eV = am3.get_egrid_lep()
    Ee_erg = Ee_eV * eV2erg
    Eg_eV = am3.get_egrid_photons()
    Eg_erg = Eg_eV * eV2erg
    Ep_eV = am3.get_egrid_had()
    Ep_erg = Ep_eV * eV2erg
    En_eV = am3.get_egrid_neutrinos()
    En_erg = En_eV * eV2erg
    return Ee_eV, Ee_erg, Eg_eV, Eg_erg, Ep_eV, Ep_erg, En_eV, En_erg