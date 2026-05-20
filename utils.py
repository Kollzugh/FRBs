# Designed to stock the argument parameters and all the functions used in our differents code.
# This file is build with the following section: arguments, plot functions, fit function, computing function and miscellaneous.

import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import corner

from scipy.ndimage import gaussian_filter
from pygdsm import GlobalSkyModel
import healpy as hp

from scipy.integrate import quad

from astropy.coordinates import SkyCoord    # convert in celestial coordinates
from astropy.table import Table             # to use vstack
import astropy.units as u                   # convert unit degrees/radians mainly
from astroquery.sdss import SDSS            # Sloan Digital Sky Survey (SDSS) catalog to cross-match

#from matplotlib.font_manager import fontManager, FontProperties
#from dl import queryClient as qc, authClient as ac

from importlib.metadata import version      # DESI software
desitarget_version = version("desitarget")


from constants import c, G, m_p, H0, Omega_m, Omega_L, Omega_b, f_IGM, Y_H, Y_He, X_e_H, X_e_He, tun_par_SDSS, tun_par_DESI, Nbr_sorted, FRB_table_file, DM_host
from data import data


#-----------ARGUMENTS---------------------------------------------------------------------------------------------
def get_args() -> argparse.Namespace:                       # load the arguments needed for the scripts
    parser = argparse.ArgumentParser()                      # first we have to create a parser object
    parser.add_argument(
        "--model",
        type=str,
        default="MP18",
        help="MP18 for Macquart relation with Planck18 values or MP18_var to vary ±10% H0, Omega_b,L",
    )
    parser.add_argument(
        "--type",
        type=str,
        default="obs",
        help="cosmic for DM_cosmic alone or obs for DM_obs = DM_Cosmic + DM_host * (1+z)^-1",
    )
    parser.add_argument(
        "--plot",
        type=str,
        default=None,
        help="mollweide or DM_histogram",
    )
    parser.add_argument(
        "--coord",
        type=str,
        default="galactic",
        help="equatorial or galactic",
    )
    parser.add_argument(
        "--catalog",
        type=str,
        default="no",
        help="To create a catalog of FRB/host catalog: yes or no",
    )

    args = parser.parse_args()
    return args
args_utils = get_args()                                            # to use the arguments in this file


#-----------PLOT FUNCTIONS----------------------------------------------------------------------------------------
def plot_DM(z, DM, fit_line, slope, intercept, DM_host_val):                              # plot theoretical DM curves
    if args_utils.type == "cosmic":
        plt.plot(z, DM, label="MP18_cosmic")
        if fit_line is not None and slope is not None and intercept is not None:
            plt.plot(z, fit_line, '--', label=fr"$Fit:\ a={slope:.2f},\ b={intercept:.2f}$")
    elif args_utils.type == "obs":
        plt.plot(z, DM, label=f"$MP18,\ DM_{{host}}={DM_host_val}pc/cm^3$" if DM_host_val is not None else "Error")
        if slope is not None and intercept is not None:
            plt.plot([], [], '--', label=fr"$Fit:\ a={slope:.2f},\ b={intercept:.2f}$")   # the fit lines are not plotted but we keep the pars values

def plot_DM_var(z, DM):                                                                   # plot DM curves with ±10% of variation on 3 pars
    plt.plot(z, DM[0], label="std")    
    plt.plot(z, DM[1], label="H0 +10%")
    plt.plot(z, DM[2], label="H0 -10%")
    plt.plot(z, DM[3], label="Ω_b +10%")
    plt.plot(z, DM[4], label="Ω_b -10%")
    plt.plot(z, DM[5], label="Ω_m +10%")
    plt.plot(z, DM[6], label="Ω_m -10%")

def mollweide(RADEC_CHIME, RADEC_PASA, RADEC_ASKAP, RADEC_REP_CHIME):                     # plot mollweide projections
    plt.figure(figsize=(10,6))
    ax = plt.subplot(projection="mollweide")

    if args_utils.coord == "galactic":                                                    # better to use galactic coord for galactic objects (if measurment values in degrees)
        gsm = GlobalSkyModel(freq_unit='MHz')
        sky_map = gsm.generate(408)
        lon = np.linspace(-np.pi, np.pi, 1000)
        lat = np.linspace(-np.pi/2, np.pi/2, 500)
        LON, LAT = np.meshgrid(lon, lat)
        pix = hp.ang2pix(
            hp.get_nside(sky_map),
            np.pi/2 - LAT,
            (-LON) % (2*np.pi)
        )
        Z_smooth = gaussian_filter(np.log10(sky_map)[pix], sigma=3)
        level = np.percentile(Z_smooth, 80)
        fig_tmp, ax_tmp = plt.subplots()
        cs = ax_tmp.contour(LON, LAT, Z_smooth, levels=[level])
        plt.close(fig_tmp)
        for path in cs.collections[0].get_paths():
            x, y = path.vertices.T
            jump = np.sqrt(np.diff(x)**2 + np.diff(y)**2) > 0.15
            cuts = np.where(jump)[0] + 1
            for xx, yy in zip(np.split(x, cuts), np.split(y, cuts)):
                if len(xx) > 20:
                    ax.plot(xx, yy, color="gray", linewidth=1.2)

        datasets_gal = [RADEC_CHIME, RADEC_PASA, RADEC_ASKAP, RADEC_REP_CHIME]
        (l_CHIME, b_CHIME), (lon_PASA, lat_PASA), (lon_ASKAP, lat_ASKAP), (lon_rep_C, lat_rep_C) = map(mollweide_galactic, datasets_gal)

        ax.scatter(l_CHIME, b_CHIME, s=10, color='red', marker='x', label="CHIME") 
        ax.scatter(lon_PASA, lat_PASA, s=10, color='dodgerblue', marker='x', label="FRBCAT") 
        ax.scatter(lon_ASKAP, lat_ASKAP, s=10, color='black', marker='x', label="ASKAP")
        ax.scatter(lon_rep_C, lat_rep_C, s=20, color='black', marker='D', label="Repeated FRB")
        ax.set_title("FRBs sky distribution (galactic coord: l-b)", pad=20)

    elif args_utils.coord == "equatorial":                                                # better to use equatorial coord for extragalactic objects
        datasets_eq = [RADEC_CHIME, RADEC_PASA, RADEC_ASKAP, RADEC_REP_CHIME]
        (RA_CHIME, DEC_CHIME), (RA_PASA, DEC_PASA), (RA_ASKAP, DEC_ASKAP), (RA_REP_C, DEC_REP_C) = [
            np.deg2rad(np.array(d).T) for d in datasets_eq
        ]
        RA_CHIME, RA_PASA, RA_ASKAP, RA_REP_C = mollweide_equtorial(RA_CHIME), mollweide_equtorial(RA_PASA), mollweide_equtorial(RA_ASKAP), mollweide_equtorial(RA_REP_C)
        ax.scatter(RA_CHIME, DEC_CHIME, s=10, color='red', marker='x', label="CHIME")
        ax.scatter(RA_PASA, DEC_PASA, s=10, color='dodgerblue', marker='x', label="FRBCAT")
        ax.scatter(RA_ASKAP, DEC_ASKAP, s=10, color='black', marker='x', label="ASKAP-")
        ax.scatter(RA_REP_C, DEC_REP_C, s=20, color='black', marker='D', label="Repeated FRB") 
        ax.set_title("FRBs sky distribution (equatorial coord: $\\alpha-\\delta$)", pad=20)
    else:
        print("You have to choose the kind of coordinates you want for the Mollweide plot (equatorial or galactic).")

    ax.grid(True)
    ax.legend(loc="lower right", bbox_to_anchor=(1.11, -0.13))

def histo_DM(DM_nonrep, DM_rep, bins=None):                                         # plot a FRBs histogram in function of DM
    if bins is None:                                                                # to have equal width for repeaters and non repeaters
        all_dm = np.concatenate([DM_nonrep, DM_rep])
        bins = np.linspace(min(all_dm), max(all_dm), 30)
    plt.hist(DM_nonrep, bins, alpha=0.5, label="Non-repeaters")
    plt.hist(DM_rep, bins, alpha=0.5, label="Repeaters")
    #plt.yscale("log")
    plt.xlabel("DM (pc/cm$^{3}$)")
    plt.ylabel("Counts")
    plt.legend(loc="upper right")

def corner_FRB(H0_samples, Omega_b_samples, f_IGM_samples):                         # Not useful because H0 is not a MCMC cloud
    corner_datasets = np.column_stack([H0_samples, Omega_b_samples, f_IGM_samples])
    labels = [r"$H_0$", r"$\Omega_b$", r"$f_{\rm IGM}$"]
    fig_corner = corner.corner(
        corner_datasets,
        labels=labels,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_fmt=".3f",
        label_kwargs={"fontsize": 14},
        title_kwargs={"fontsize": 12},
        #plot_datapoints=True,
        #fill_contours=False,
        levels=(0.68, 0.95),
        truths=[
            67.66,    # H0 Planck
            0.0488,    # Omega_b Planck environ
            0.84,     # f_IGM
        ]
    )
    fig_corner.tight_layout()
    plt.show()

#-----------FIT FUNCTIONS-----------------------------------------------------------------------------------------
def linear_fit(x, y):                                                           # compute a, b and fit for a line: y = a.x + b
    p, cov = np.polyfit(x, y, 1, cov=True)
    slope, intercept = p
    dslope, dintercept = np.sqrt(np.diag(cov))
    fit_line = slope * x + intercept
    return fit_line, slope, intercept, dslope, dintercept

def log_quadratic_fit(z, DM):
    z, DM = np.asarray(z), np.asarray(DM)
    mask = (z > 0) & (DM > 0)
    z_fit_data, DM_fit_data = z[mask], DM[mask]
    a, b, c = np.polyfit(np.log10(z_fit_data), np.log10(DM_fit_data), 2)
    return a, b, c


#-----------COMPUTING FUNCTIONS-----------------------------------------------------------------------------------
def E(z, Omega_m_val):                                                          # to compute E in the Macquart relation (see "Project_intro.pdf")
    E_val = (Y_H * X_e_H + 0.5 * Y_He * X_e_He) / np.sqrt(Omega_m_val * (1+z)**3 + Omega_L)
    return E_val

def integral(z, Omega_m_val):                                                   # define the integral of the Macquart relation
        integrand = lambda z: (1+z) * E(z, Omega_m_val)
        val, _ = quad(integrand, 0, z)
        return val

def DM_cosmic(z_FRB, model):                                                    # compute DM with the Macquart relation
    factor_std = (f_IGM * Omega_b * (3 * H0 * c)) / (8 * np.pi * G * m_p)
    if model == "MP18":                                                         # plot only the line with Planck 2018 coef
        I = integral(z_FRB, Omega_m)
        return factor_std * I
    elif model == "MP18_var":                                                   # plot lines with Planck 2018 coef and the variations of ±10%
        I_std   = integral(z_FRB, Omega_m)
        I_Omp10 = integral(z_FRB, Omega_m * 1.1)
        I_Omm10 = integral(z_FRB, Omega_m * 0.9)
        fac_H0p10 = (f_IGM * Omega_b * (3 * (H0*1.1) * c)) / (8 * np.pi * G * m_p)
        fac_H0m10 = (f_IGM * Omega_b * (3 * (H0*0.9) * c)) / (8 * np.pi * G * m_p)
        fac_Obp10 = (f_IGM * (Omega_b*1.1) * (3 * H0 * c)) / (8 * np.pi * G * m_p)
        fac_Obm10 = (f_IGM * (Omega_b*0.9) * (3 * H0 * c)) / (8 * np.pi * G * m_p)
        return (
            factor_std * I_std,      # standard Planck 2018
            fac_H0p10 * I_std,       # H0 +10%
            fac_H0m10 * I_std,       # H0 -10%
            fac_Obp10 * I_std,       # Omega_b +10%
            fac_Obm10 * I_std,       # Omega_b -10%
            factor_std * I_Omp10,    # Omega_m +10%
            factor_std * I_Omm10     # Omega_m -10%
        )
    else:
        raise ValueError("Unknown model. You have to specify the model: MP18 for std values, or MP18_var to add ±10% variations.")

def mollweide_equtorial(ra):                                                        # useful to get equatorial coord
    ra = np.remainder(ra + 2*np.pi, 2*np.pi)                                        # force FRBs to be in [0,2pi], (modulo 2pi)
    ra[ra > np.pi] -= 2*np.pi                                                       # move to the interval [-pi,pi] if not
    return -ra                                                                      # the minus sign to invert RA (it must increase towards the left)

def mollweide_galactic(RADEC):                                                      # useful to get galactic coord
    ra_deg = [x for x, y in RADEC]
    dec_deg = [y for x, y in RADEC]
    c = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')                  # create a sky coordinates object
    l = -c.galactic.l.wrap_at(180*u.deg).radian                                     # longitude
    b = c.galactic.b.radian                                                         # latitude
    return l, b

def DM_z_real():
    combi = Table.read(FRB_table_file, format="ascii.basic", guess=False)           # read the file created from the cross-match
    unique_rows = []                                                                # this line is to avoid having the same match twice (eg, problematic for a fit)
    seen = set()
    for row in combi:
        name = row['FRB_name']
        if name not in seen:
            unique_rows.append(row)
            seen.add(name)
    combi = Table(rows=unique_rows, names=combi.colnames)

    DM = np.array(combi['DM_obs'], dtype=float)
    DM_err = np.array(combi['DM_obs_err'], dtype=float)
    z = np.array(combi['z'], dtype=float)
    z_err = np.array(combi['z_err'], dtype=float)
    
    a, b = np.polyfit(z, DM, 1)
    z_fit = np.linspace(min(z), max(z), 200)
    DM_fit = a*z_fit + b
    return z, z_err, DM, DM_err, z_fit, DM_fit, a, b


def calculate_omegab_fIGM(z, DM, ax_top):
    H0 = np.linspace(67.66, 73.00, 10) * (3.24076e-20)
    omegab, fIGM = np.array([]), np.array([])
    for i in range(len(DM_host)):
        for j in range(len(H0)):
            DM_cosmic_Y = (np.array(DM) - (DM_host[i] / (1 + z))) * 3.0857e22         # need to be in pc/m^3
            I = np.array([integral(zi, Omega_m) for zi in z])
            X = ((3 * H0[j] * c * f_IGM) / (8 * np.pi * G * m_p)) * I                 # calculation for omega_b
            fit_line, slope, intercept, slop_err, inter_err = linear_fit(X, DM_cosmic_Y)
            XX = ((3 * H0[j] * c * Omega_b) / (8 * np.pi * G * m_p)) * I              # calculation for f_IGM
            fit_lineB, slopeB, interceptB, slop_errB, inter_errB = linear_fit(XX, DM_cosmic_Y)
            omegab = np.append(omegab, slope)
            fIGM = np.append(fIGM, slopeB)

            if j == 0:
                ax_top.text(0.04, 0.93, f"$H_0$ = 67.66 km/s/Mpc: $\\Omega_{{b-}}$ = {slope:.4f}±{slop_err:.4f} and $f_{{IGM}}$ = {slopeB:.4f}±{slop_errB:.4f}", transform=ax_top.transAxes, fontsize=10, va='top')
            elif j == 9:
                ax_top.text(0.04, 0.86, f"$H_0$ = 73.00 km/s/Mpc: $\\Omega_{{b+}}$ = {slope:.4f}±{slop_err:.4f} and $f_{{IGM}}$ = {slopeB:.4f}±{slop_errB:.4f}", transform=ax_top.transAxes, fontsize=10, va='top')
    return H0/3.24076e-20, DM_cosmic_Y, omegab, fIGM

def variance(z, res_DM):
    bins = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0])                                   # to compute res vs variance
    bin_id = np.digitize(z, bins)
    z_bin_center, var_bin = [], []
    for k in range(1, len(bins)):
        mask = bin_id == k
        if np.sum(mask) > 1:
            z_bin_center.append(0.5 * (bins[k-1] + bins[k]))
            var_bin.append(np.var(res_DM[mask], ddof=1))
    plt.plot(z_bin_center, var_bin, ".-")
    plt.xlabel("z")
    plt.ylabel(r"Variance")
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.show()

def create_FRB_to_galaxie_catalog():                                                # to create your own FRB/galaxie matches catalog (it can take a while ~5-10mn)
    coords_all, names, dm, dmE,raE, decE = merge_and_skycoord(                      # if you prefer, the catalog was already created in data/FRBs_galaxies_catalog.txt
        data.radec_CHIME, data.radec_PASA, data.radec_ASKAP, 
        data.names_CHIME, data.names_PASA, data.names_ASKAP, 
        data.dm_CHIME, data.dm_PASA,
        data.raErr_CHIME, data.decErr_CHIME,
        data.dm_err_CHIME, data.dm_err_PASA,
        )

    CrossMatch_SDSS = cross_match_SDSS(coords_all, tun_par_SDSS)
    CrossMatch_DESI = cross_match_DESI(coords_all, tun_par_DESI, Nbr_sorted)
    matches = []
    for i in range(len(coords_all)):

        r_sdss = CrossMatch_SDSS[i]
        if r_sdss is not None and len(r_sdss) > 0:
            row = {}
            row['FRB_name'] = names[i]
            row['ra'] = r_sdss['ra'][0]
            row['ra_err'] = str(raE[i])
            row['dec'] = r_sdss['dec'][0]
            row['dec_err'] = str(decE[i])
            row['DM_obs'] = str(dm[i])
            row['DM_obs_err'] = str(dmE[i])
            row['z'] = float(r_sdss['z'][0]) if 'z' in r_sdss.colnames else None
            row['z_err'] = float(r_sdss['zErr'][0]) if 'zErr' in r_sdss.colnames else None
            row['class'] = str(r_sdss['class'][0]) if 'class' in r_sdss.colnames else None
            row['subClass'] = str(r_sdss['subClass'][0]) if 'subClass' in r_sdss.colnames else None
            row['catalog'] = 'SDSS'
            if row['z'] is not None and row['z'] > 0.01:
                matches.append(row)

        r_desi = CrossMatch_DESI[i]
        if r_desi is not None and len(r_desi) > 0:
            row = {}
            row['FRB_name'] = names[i]
            row['ra'] = r_desi['ra'].iloc[0]
            row['ra_err'] = str(raE[i])
            row['dec'] = r_desi['dec'].iloc[0]
            row['dec_err'] = str(decE[i])
            row['DM_obs'] = str(dm[i])
            row['DM_obs_err'] = str(dmE[i])
            row['z'] = float(r_desi['z'].iloc[0]) if 'z' in r_desi.columns else None
            row['z_err'] = float(r_desi['zerr'].iloc[0]) if 'zerr' in r_desi.columns else None
            row['class'] = str(r_desi['spectype'].iloc[0]) if 'spectype' in r_desi.columns else None
            val = r_desi['subtype'].iloc[0] if 'subtype' in r_desi.columns else None
            row['subClass'] = '--' if pd.isna(val) else str(val)            
            row['catalog'] = 'DESI'
            if row['z'] is not None and row['z'] > 0.01:
                matches.append(row)

    if len(matches) > 0:
        combined = Table(rows=matches)
        combined = combined['FRB_name', 'ra', 'ra_err', 'dec', 'dec_err',
                    'DM_obs', 'DM_obs_err', 'z', 'z_err',
                    'class', 'subClass', 'catalog']
        combined.pprint(max_lines=-1, max_width=-1)
        combined.write("data/FRBs_galaxies_catalog.txt", format="ascii", overwrite=True)
        return combined

#-----------MISCELLANEOUS-----------------------------------------------------------------------------------------
def merge_and_skycoord(RADEC_CHIME, RADEC_PASA, RADEC_ASKAP, NAMES_CHIME, NAMES_PASA, NAMES_ASKAP, DM_CHIME, DM_PASA, RA_ERR, DEC_ERR, DM_ERR_CHIME, DM_ERR_PASA):
    radec_all = list(RADEC_CHIME) + list(RADEC_PASA) + list(RADEC_ASKAP)
    names = list(NAMES_CHIME) + list(NAMES_PASA) + list(NAMES_ASKAP)
    dm = list(DM_CHIME) + list(DM_PASA) + [np.nan]*len(RADEC_ASKAP)                 # note that no data of DM is available for ASKAP
    dm_err = list(DM_ERR_CHIME) + list(DM_ERR_PASA) + [np.nan]*len(RADEC_ASKAP)
    ra_err = list(RA_ERR) + [np.nan]*len(RADEC_PASA) + [np.nan]*len(RADEC_ASKAP)
    dec_err = list(DEC_ERR) + [np.nan]*len(RADEC_PASA) + [np.nan]*len(RADEC_ASKAP)
    ra, dec = zip(*radec_all)
    coordinates = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    return  coordinates, names, dm, dm_err, ra_err, dec_err                         # create an object for a cross-match and keep FRB names

def cross_match_SDSS(COORDS_ALL, tun_SDSS):                                         # to take data from the SDSS catalog
    results = []
    for coord in COORDS_ALL:
        res = SDSS.query_region(                                                    # specific query
            coord,
            radius=tun_SDSS*u.arcsec,
            spectro=True,
            photoobj_fields=['ra', 'dec'],
            specobj_fields=['z', 'zErr', 'class', 'subClass']
        )
        results.append(res)
    return results

def cross_match_DESI(COORDS_ALL, tun_DESI, top_n):                                  # to cross-match with the DESI catalog
    results = []
    radius_deg = (tun_DESI * u.arcsec).to(u.deg).value
    for coord in COORDS_ALL:
        ra = coord.ra.deg
        dec = coord.dec.deg
        query = f"""
        SELECT TOP {top_n}
            p.targetid,
            p.ra, p.dec,
            z.z, z.zerr, z.spectype, z.subtype
        FROM desi_edr.photometry AS p
        JOIN desi_edr.zpix AS z
            ON p.targetid = z.targetid
        WHERE q3c_radial_query(p.ra, p.dec, {ra}, {dec}, {radius_deg})
          AND z.z IS NOT NULL
        """
        try:
            res = qc.query(sql=query, fmt='pandas', timeout=600)
        except Exception as e:
            print(f"Erreur pour RA={ra}, DEC={dec}: {e}")
            res = pd.DataFrame()
        results.append(res)
    return results