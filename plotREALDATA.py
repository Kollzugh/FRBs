# Designed play with real data.

import numpy as np
import matplotlib.pyplot as plt

from constants import DM_host, DM_host_err, DM_MW
from data import data
from utils import get_args, histo_DM, mollweide, DM_z_real, linear_fit, create_FRB_to_galaxie_catalog, variance, log_quadratic_fit, calculate_omegab_fIGM, corner_FRB
args = get_args()


#-----------MAIN CODE---------------------------------------------------------------------------------------------
if args.plot == "mollweide":
    mollweide(data.radec_CHIME, data.radec_PASA, data.radec_ASKAP, data.radec_rep_CHIME)
elif args.plot == "DM_histogram":
    dm_nonrep = np.concatenate((data.dm_PASA, data.dm_nonrep_CHIME))                    # dm_ASKAP not furnished
    dm_rep = data.dm_rep_CHIME                                                          # rep from FRBCAT/PASA already in CHIME
    histo_DM(dm_nonrep, dm_rep)
elif args.plot == "DM_z":
    z_vals, z_vals_err, DM_obs, DM_obs_err, z_fit, DM_fit, a, b = DM_z_real()
    if args.catalog == "yes":                                                           # create a catalog of candidats after a cross-match
        create_FRB_to_galaxie_catalog()                                                 # change tunable parameters to select the area in constants.py

    fig = plt.figure(figsize=(7, 7))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 1.2], hspace=0.05)
    ax_top = fig.add_subplot(gs[0, 0])
    ax_bot = fig.add_subplot(gs[1, 0], sharex=ax_top)

    H0, DM_IGM, calc_omegab, calc_fIGM = calculate_omegab_fIGM(z_vals, DM_obs, ax_top)  # calculate with the Macquart relation
    #corner_FRB(H0, calc_omegab, calc_fIGM)                                             # plot to find best values of H0, Omega_b and f_IGM > not useful (see the function in utils.py)

    """#UNCOMENT this part TO PLOT omegab/fIGM vs H0
    fig_tmp, ax_tmp = plt.subplots(1, 2, figsize=(10, 4))    
    ax_tmp[0].plot(H0, calc_omegab, marker='.', ls=' ')
    ax_tmp[1].plot(H0, calc_fIGM, marker='.', ls=' ')
    ax_tmp[0].set_ylabel("$\\Omega_{{b}}$")
    ax_tmp[1].set_ylabel("$f_{{IGM}}$")
    ax_tmp[0].set_xlabel("$H_0$ (km/s/Mpc)")
    ax_tmp[1].set_xlabel("$H_0$ (km/s/Mpc)")
    fig_tmp.show()"""

    ############# Fit #############
    fit_line, slope_lin, intercept_lin, slop_err_lin, inter_err_lin = linear_fit(z_vals, DM_IGM/3.0857e22)  # linear model
    z_line = np.linspace(0, 1.2, 300)
    DM_line = slope_lin * z_line + intercept_lin 

    A, B, C = log_quadratic_fit(z_vals, DM_IGM/3.0857e22)                                                   # quadratic model
    z_quad = np.linspace(1e-4, 1.2, 200)
    DM_quad = 10**(A * (np.log10(z_quad))**2 + B * np.log10(z_quad) + C)

    ############# Uncertainty Bands #############
    band = 100 + z_line*250                                                                                 # linear model
    DM_upper = DM_line + band
    DM_lower = DM_line - band                                                        

    band_quad = 10**(-0.15*z_quad**2 + 0.8*z_quad + 2)                                                      # quadratic model
    DM_upper_quad = DM_quad + band_quad
    DM_lower_quad = DM_quad - band_quad

    ############# Residuals and metrics #############
    DM_model_lin = slope_lin * z_vals + intercept_lin                                                       # residuals on DM_IGM for the linear model
    residuals_lin = DM_IGM/3.0857e22 - DM_model_lin
    residuals_err_lin = np.sqrt(DM_obs_err**2 + (DM_host_err / (1 + z_vals))**2)
    DM_model_quad = 10**(A * (np.log10(z_vals))**2 + B * np.log10(z_vals) + C)                              # from the quadratic model but useful to calculate the metrics

    x_lin = np.median(np.abs(DM_obs - (DM_model_lin + DM_MW + DM_host[0]/(1+z_vals))))                      # calculation of the median error x
    print('Linear model: x_lin = ', x_lin)
    x_quad = np.median(np.abs(DM_obs - (DM_model_quad + DM_MW + DM_host[0]/(1+z_vals)))) 
    print('Quadratic log model: x_quad = ', x_quad)
    RMSE_lin = np.sqrt((1/len(z_vals)) * np.sum((DM_obs - (DM_model_lin + DM_MW + DM_host[0]/(1+z_vals)))**2)) # calculation of the RMSE
    print('Linear model: RMSE_lin = ', RMSE_lin)
    RMSE_quad = np.sqrt((1/len(z_vals)) * np.sum((DM_obs - (DM_model_quad + DM_MW + DM_host[0]/(1+z_vals)))**2))
    print('Quadratic log model: RMSE_quad = ', RMSE_quad)

    ############# Top graph: DM(z #############
    ax_top.errorbar(z_vals, DM_obs, xerr=z_vals_err, yerr=DM_obs_err, fmt='.', ms=4, capsize=3, ecolor='black')
    ax_top.fill_between(z_line, DM_lower, DM_upper, alpha=0.15, color='orange')
    ax_top.plot(z_line, DM_line, lw=2, label=f"DM = {slope_lin:.1f} z + {intercept_lin:.1f}")
    ax_top.fill_between(z_quad, DM_lower_quad, DM_upper_quad, alpha=0.15, color='green')
    ax_top.plot(z_quad, DM_quad, lw=2, color='green', label=fr"$\log DM = {A:.2f}(\log z)^2 + {B:.2f}\log z + {C:.2f}$")
    ax_top.set_xlim(0, 1.2)
    ax_top.set_ylim(0, 2000)
    ax_top.set_ylabel("DM$_{IGM}$ (pc cm$^{-3}$)")
    ax_top.grid(True, which="both", alpha=0.7)
    ax_top.legend(loc="upper left", bbox_to_anchor=(0.01, 0.8), fontsize=10)
    
    ############# Bottom graph: Res(z) #############
    ax_bot.errorbar(z_vals, residuals_lin, yerr=residuals_err_lin, fmt='.', ms=4, capsize=3, ecolor='black', label="Residuals")
    ax_bot.axhline(0, ls='--', color='gray')
    ax_bot.set_ylabel("Residuals (DM$_{IGM}$)")
    ax_bot.set_xlabel("z")
    ax_bot.grid(True, which="both", alpha=0.7)
    plt.setp(ax_top.get_xticklabels(), visible=False)

plt.show()

"""#UNCOMENT the following lines TO PLOT the VARIANCE and calculate chi square:
res_obs = DM_obs - ((slope_lin * z_vals + intercept_lin) + (DM_host/(1+z_vals)) + DM_MW)    # to plot the variance vs z
chi_square = (1 / (len(z_vals)-2)) * np.sum((res_obs / DM_obs_err)**2)                      # DM_obs! huge value: normal (not used for the variance)
variance(z_vals, res_obs)"""