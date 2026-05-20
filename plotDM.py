# Designed to plot DM.

import numpy as np
import matplotlib.pyplot as plt

from utils import get_args, DM_cosmic, linear_fit, plot_DM, plot_DM_var
args = get_args()


#-----------MAIN CODE---------------------------------------------------------------------------------------------
DM_host = [0, 25, 50, 100]                                                          # pc cm^-3
z_vals = np.linspace(0, 1.3, 100)                                                   # simulate values of z to calculate DM

DM_Macquart = [DM_cosmic(z, model=args.model) for z in z_vals]                      # compute DM from the Macquart relation
DM_Macquart_conv = np.array(DM_Macquart) / 3.0857e22                                # put DM in pc cm^-3

if args.model == "MP18":                                                            # plot the Macquart relation with Planck18 values
    if args.type == "cosmic":                                                       # if only DM_IGM is needed
        fit_line_MP18, slope_MP18, intercept_MP18, _, _ = linear_fit(z_vals, DM_Macquart_conv)    
        plot_DM(z_vals, DM_Macquart_conv, fit_line_MP18, slope_MP18, intercept_MP18, None)
        plt.ylabel("DM$_{Cosmic}$ (pc cm$^{-3}$)")
    elif args.type == "obs":                                                        # if we need DM_obs = DM_cosmic(z) + DM_host(z) (without DM_MW !)
        DM_OBS = np.zeros(len(DM_host))
        for i in range(len(DM_host)):
            DM_OBS = DM_Macquart_conv + (DM_host[i] / (1 + z_vals))
            fit_line_MP18_obs, slope_MP18_obs, intercept_MP18_obs, _, _ = linear_fit(z_vals, DM_OBS)    
            plot_DM(z_vals, DM_OBS, fit_line_MP18_obs, slope_MP18_obs, intercept_MP18_obs, DM_host[i])
        plt.ylabel("DM$_{obs}$ (pc cm$^{-3}$)")
    else:
        raise ValueError("Unknown type. You have to specify the type of relation: cosmic for DM_cosmic alons, or obs for DM_obs = DM_Cosmic + DM_host * (1+z)^-1.")
elif args.model == "MP18_var":                                                      # compute DM from the Macquart relation with variable parameters
    DM_Macquart_conv = DM_Macquart_conv.T
    plot_DM_var(z_vals, DM_Macquart_conv)
    plt.ylabel("DM$_{Cosmic}$ (pc cm$^{-3}$)")

plt.legend(loc="upper left", fontsize=9)
plt.xlim(0, 1.2); plt.ylim(0, 1200)
plt.xlabel("z")
plt.grid()
plt.show()
