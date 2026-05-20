# Designed to stock the constants used in our project.

#------------GENERAL------------
c = 2.9979e8                                        # m/s
G = 6.6738e-11                                      # m^3/(kg s^2)
m_p = 1.6726e-27                                    # kg

#------------DM_Macquart------------
H0 = 67.66 * (3.24076e-20)                          # 1/s because we take all other values in SI
Omega_m = 0.3103                                    # matter density parameter (total matter: dark matter + baryons), Planck18
Omega_L = 0.6897                                    # dark energy density parameter (cosmological constant), Planck18
Omega_b = 0.0488                                    # baryon density parameter (ordinary matter), Planck18

f_IGM = 0.84                                        # fraction of baryonic matter in the intergalactic medium (IGM)
Y_H = 0.75                                          # mass fraction of hydrogen
Y_He = 0.25                                         # mass fraction of helium
X_e_H = 1.0                                         # ionization fraction of hydrogen (1 = fully ionized)
X_e_He = 1.0                                        # ionization fraction of helium (1 = fully ionized)

#------------General DM-------------
DM_host = [50]                                      # pc cm^-3, you can put other values to compare
DM_host_err = 50
DM_MW = 80

#------------Cross-match------------
tun_par_SDSS = 20                                   # useful to tune the crossmatch
tun_par_DESI = 45
Nbr_sorted = 20                                     # selection of the max number of target to sort out of the DESI catalog
FRB_table_file = "data/FRBs_candidats_catalog.txt"  # to use ASKAP data