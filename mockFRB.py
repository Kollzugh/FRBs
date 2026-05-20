# Designed to create a mock FRB catalog.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#-----------VARIABLES---------------------------------------------------------------------------------------------

Nbr_FRB = 100

z_mean = 0.75
z_std_dev = 0.3

#-----------MAIN CODE---------------------------------------------------------------------------------------------
np.random.seed(42)
z = np.random.normal(z_mean, z_std_dev, Nbr_FRB)
z = np.clip(z, 0.01, None)                          # To avoid z=0 or z<0

DM_IGM = 1000 * z                                   # To add a weight on DM_IGM
cosmic_noise = np.random.normal(0, 200, Nbr_FRB)    # To simulate the cosmic web fluctuations: dispersion ~100 pc/cm^3
DM_host = np.random.normal(100, 30, Nbr_FRB) / (1 + z)
#DM_MW = np.random.normal(50, 10, Nbr_FRB)          # Not added because not included in question 6 of "Theoretical framework and modeling"

DM_total = DM_IGM + cosmic_noise + DM_host          # without DM_MW, add it if nexessary

catalog = pd.DataFrame({
    "Redshift": z,
    "DM_IGM": DM_IGM,
    "DM_obs": DM_total
})

#print(catalog.head())
plt.scatter(z, DM_total, s=10, alpha=0.6)
plt.xlim(0, 1.2)
plt.ylim(0, 1200)
plt.xlabel("z")
plt.ylabel("$DM_{{obs}} (pc/cm^{{3}})$")
plt.title("Mock FRB Catalog")
plt.show()