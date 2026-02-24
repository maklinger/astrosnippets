import numpy as np
from .colorschemes import tol_cset


def Eobs(E_co, z=0, doppler=1):
    """
    converts comoving to observed energy with average angle cos(theta) = beta
    args:
        E_co   : comoving energy [any energy unit]
        z      : redshift
        doppler: Doppler factor of blob
    returns:
        observed energy [same unit as E_co]
    """
    return E_co * doppler / (1+z)


def dens2flux(E2dn_dE, doppler, tesc, V, distance, doppler_power=4):
    """
    converts comoving density into flux [erg/s/cm^2] inclduing doppler factor
    args:
        E2dn_dE : comov. energy spectrum of free streaming particles [erg]
        doppler: Doppler factor of blob
        t_esc   : dynamic time scale [s]
        V       : volume of blob [cm³]
        distance: distance to blob [cm]
        doppler_power: exponent of dopplerfactor, defaults to blob value of 4
    returns:
        observed energy flux [erg/(cm^2s)]
    """
    return doppler**doppler_power * E2dn_dE/tesc * V /(4*np.pi*distance**2)


# for absorption ranges with fading
def plot_logfaderange(ax, E1, E2, facE1, facE2, color, N=5):
    # main range alpha=1
    ax.axvspan(E1, E2, color=color, alpha=1)
    # N bins to fade out
    bounds1 = np.logspace(np.log10(E1/facE1), np.log10(E1), N)
    bounds2 = np.logspace(np.log10(E2), np.log10(E2*facE2), N)
    for i in range(N-1):
        ax.axvspan(bounds1[i], bounds1[i+1], color=color, alpha=(i+1)/N)
        ax.axvspan(bounds2[i], bounds2[i+1], color=color, alpha=1-(i+1)/N)
    return ax


def plot_spectrum(
    ax, doppler, distance, z, t_esc, V, Eg_eV, En_eV,
    E2dN_dE_tot_gamma, E2dN_dE_tot_nu, components,
    # grey bands for absorption
    cgrey="#444444", Emins_grey=[1, 1e13], Emaxs_grey=[3e2, 1e22], facE_grey=5, Ngrey=30,
    labels_grey=["dust\nphotoel.", r"$\gamma \gamma$ (EBL,CMB)"], textcolgrey="lightgrey",
    # boxes indicating observational ranges
    Emins_obs=[1e-4, 1e3, 3e11], Emaxs_obs=[1e-3, 1e5, 1e13],
    dlgobs=-0.7, color_obs="lightgrey", labels_obs=["radio", "X-ray", "VHE"],
    alpha_dom=0.4, zorder_dom=10, lw_main=3, col_tot="k", col_nu="teal",
    ylim=[1e-18, 1e-5], xlim=[1e-4, 1e18], doppler_power=4, lines={}
):

    # grey band for absorption
    ygrey = 0.8 * ylim[1]
    for i in range(len(Emins_grey)):
        ax = plot_logfaderange(
            ax, Emins_grey[i], Emaxs_grey[i], facE_grey, facE_grey, cgrey, Ngrey)
        # label for grey band
        ax.text(np.sqrt(Emins_grey[i]*Emaxs_grey[i]), ygrey, labels_grey[i],
                va="top", ha="center", ma="center", color=textcolgrey)

    # observational ranges
    ymax_obs = ylim[1]
    ymin_obs = 10**dlgobs * ymax_obs
    for i in range(len(Emins_obs)):
        ax.fill_between([Emins_obs[i], Emaxs_obs[i]], ymin_obs*np.ones(2),
                        ymax_obs*np.ones(2), alpha=1, color=color_obs)
        ax.text(np.sqrt(Emins_obs[i]*Emaxs_obs[i]), np.sqrt(ymin_obs*ymax_obs),
                labels_obs[i], va="center", ha="center", ma="center", color="k")

    Egobs = Eobs(Eg_eV, z, doppler)

    # total photons
    lines["total_photons"] = ax.loglog(
        Egobs, dens2flux(E2dN_dE_tot_gamma, doppler, t_esc, V, distance, doppler_power), lw=lw_main,
        zorder=zorder_dom, color=col_tot)
    
    # total neutrinos
    lines["total_neutrinos"] = ax.loglog(
        Eobs(En_eV, z, doppler),
        dens2flux(E2dN_dE_tot_nu, doppler, t_esc, V, distance, doppler_power),
        lw=lw_main, zorder=zorder_dom, color=col_nu)

    # components

    for E2dN_dE, col, dom, zup, ls, lbl, lw in components:
        if dom:
            ax.fill_between(
                Egobs, np.zeros_like(Egobs),
                dens2flux(E2dN_dE, doppler, t_esc, V, distance, doppler_power),
                color=col, alpha=alpha_dom, zorder=zorder_dom-1+zup, ls=ls,
                lw=lw)
        lines[lbl] = ax.loglog(
            Egobs, dens2flux(E2dN_dE, doppler, t_esc, V, distance, doppler_power),
            color=col, ls=ls, zorder=zorder_dom+zup, label=lbl, lw=lw)

    # title_string = r"$t=" + f"{tobs}"
    # title_string += r"$ s, $z=%g" % (z)
    # if Ekiniso is not None:
    #     title_string += r",\; E_\mathrm{kin,iso}=%.0e \; \mathrm{erg}" %(Ekiniso)
    # title_string += "$"
    # ax.set_title(title_string, loc="left")
    ax.set_ylim(*ylim)
    ax.set_yticks(10**np.arange(np.log10(ax.get_ylim()
                  [0]), np.log10(ax.get_ylim()[1])+0.5))
    ax.set_xlim(*xlim)
    ax.set_xticks(10**np.arange(-6, np.log10(ax.get_xlim()[1])+0.5, 3))
    ax.set_xticks(10**np.arange(np.log10(ax.get_xlim()
                  [0]), np.log10(ax.get_xlim()[1])+0.5, 1), minor=True)
    ax.set_xticklabels(["" for i in np.arange(
        np.log10(ax.get_xlim()[0]), np.log10(ax.get_xlim()[1])+0.5, 1)], minor=True)
    ax.set_aspect("equal")
    ax.grid(which="both", alpha=0.3, zorder=-10)
    ax.set_axisbelow(True)
    ax.set_xlabel(r"$E$ [eV]")
    ax.set_ylabel(r"$EF_E$ [erg/(cm²s)]")
    return ax, lines



def plot_coolingtimes(
        ax, Ep_eV, Ee_eV, Eg_eV, ind, Eemin_eV, Epmin_eV,
        tpsys, tpads, tpacs, tpics, tpbhs, tppgs, tppps,
        tesys, teics, teacs, tgescs, tgpairs, tgssas, tpidec, tmudec,
        cel="tab:blue", cproton="tab:purple", cgamma="tab:green", cpion="tab:cyan",
        cmuon="tab:grey", cadi="k", alpha_inj=0.2,
        ylim=[1e-1, 1e15], xlim=[1e-4, 1e21], xmin_label=1e3, leptonic=False,
        pions=True, muons=True, lines={}):

    # injection ranges
    tels = 1/(1/tesys[ind] + 1/tpads[ind][0] + 1/teics[ind])
    Eemax_eV = Ee_eV[np.argmin((tels-teacs[ind])**2)]
    tps = 1/(1/tpsys[ind] + 1/tpads[ind][0] + 1/tpics[ind] +
             1/tpbhs[ind] + 1/tppgs[ind] + 1/tppps[ind])
    Epmax_eV = Ep_eV[np.argmin((tps-tpacs[ind])**2)]
    if not leptonic:
        lines["Epinj"] = ax.axvspan(Epmin_eV, Epmax_eV, color=cproton, alpha=alpha_inj)
        lines["Epmax"] = ax.axvline(Epmax_eV, color=cproton)
    lines["Eeinj"] = ax.axvspan(Eemin_eV, Eemax_eV, color=cel, alpha=alpha_inj)
    lines["Eemax"] = ax.axvline(Eemax_eV, color=cel)
    lines["Eemin"] = ax.axvline(Eemin_eV, color=cel)
    if not leptonic:
        lines["Epmin"] = ax.axvline(Epmin_eV, color=cproton, ls="--")

    if not leptonic:
        massp = 0.938  # proton mass / GeV
        masspipm = 0.1396
        massmu = 0.10566
        if muons:
            lines["tmusy"] = ax.loglog(Ep_eV, (massmu/massp)**4 *
                    tpsys[ind], label="mu syn", c=cmuon)
            lines["tmudec"] = ax.loglog(Ep_eV, tmudec, label="mu dec", c=cmuon, ls=":")
        if pions:
            lines["tpisy"] = ax.loglog(Ep_eV, (masspipm/massp)**4 *
                    tpsys[ind], label="pi syn", c=cpion)
            lines["tpidec"] = ax.loglog(Ep_eV, tpidec, label="pi dec", c=cpion, ls=":")

    # esc
    lines["tges"] = ax.axhline(tgescs[ind][0], label="g esc", c=cgamma)

    # adi
    if leptonic:
        cadi = cel
    lines["tpad"] = ax.axhline(tpads[ind][0], label="e adi", c=cadi, ls="-")
    # if not leptonic:
    #     ax.axhline(tpads[ind][0], label="p adi", c=cproton, ls="--")

    # acc
    lines["teac"] = ax.loglog(Ee_eV, teacs[ind], label="e acc", c=cel, ls="-")
    if not leptonic:
        lines["tpac"] = ax.loglog(Ep_eV, tpacs[ind], label="p acc", c=cproton, ls="--")

    lines["tesy"] = ax.loglog(Ee_eV, tesys[ind], label="e syn", c=cel, ls="-")
    lines["teic"] = ax.loglog(Ee_eV, teics[ind], label="e ic", c=cel, ls="--")

    if not leptonic:
        lines["tpsy"] = ax.loglog(Ep_eV, tpsys[ind], label="p syn", c=cproton)
        lines["tpic"] = ax.loglog(Ep_eV, tpics[ind], label="p ic", c=cproton, ls="--")
        lines["tpbh"] = ax.loglog(Ep_eV, tpbhs[ind], label="p bh", c=cproton, ls=":")
        lines["tppg"] = ax.loglog(Ep_eV, tppgs[ind], label="p pg", c=cproton, ls="-.")
        lines["tppp"] = ax.loglog(Ep_eV, tppps[ind], label="p pp", c=cproton, ls="-")

    lines["tgpair"] = ax.loglog(Eg_eV, tgpairs[ind], label="g pair", c=cgamma, ls="--")
    lines["tgssa"] = ax.loglog(Eg_eV, tgssas[ind], label="p ssa", c=cgamma, ls="-.")

    ax.set_ylim(*ylim)
    ax.set_yticks(10**np.arange(np.log10(ax.get_ylim()
                  [0]), np.log10(ax.get_ylim()[1])+0.5))
    ax.set_xlim(*xlim)
    ax.set_xticks(10**np.arange(np.log10(xmin_label), np.log10(ax.get_xlim()[1])+0.5, 3))
    ax.set_xticks(10**np.arange(np.log10(ax.get_xlim()
                  [0]), np.log10(ax.get_xlim()[1])+0.5, 1), minor=True)
    ax.set_xticklabels(["" for i in np.arange(
        np.log10(ax.get_xlim()[0]), np.log10(ax.get_xlim()[1])+0.5, 1)], minor=True)
    ax.set_aspect("equal")
    ax.grid(which="both", alpha=0.3, zorder=-10)
    ax.set_axisbelow(True)
    ax.set_xlabel(r"$E'$ [eV]")
    ax.set_ylabel(r"$\tau'$ [s]")
    return ax, lines



def get_comp_cols():
    comps = ["mu", "in", "bh",  "pi0", "pg","pp", "pi", "p", "pair", "nu", "total"]
    # comps = ["in", "pair", "pg", "pp", "bh", "p", "pi", "mu", "pi0", "nu", "total"]
    cols = list(tol_cset("muted"))
    compcols = {}
    for i in range(len(comps)):
        compcols[comps[i]] = cols[i]
    return compcols
