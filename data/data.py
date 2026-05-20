# Designed to load/clean data.

import pandas as pd
import warnings
import os

try:
    from cfod import catalog                # try to import the CHIME/FRB catalog (could be problematic)
except Exception as e:
    catalog = None
    print(f"cfod catalog unavailable, using local CHIME file instead. ({e})")                

from astroquery.vizier import Vizier        # to import the FRBCAT catalog from VizieR
from astropy.coordinates import SkyCoord    # to convert in celestial coordinates
import astropy.units as u                   # to convert unit degrees/radians mainly

#-----------FUNCTIONS---------------------------------------------------------------------------------------------
def extract_ra_dec(filepath):                                            # read each header file and extract data (ASKAP catalog have to be downloaded by hand)
    ra, dec = None, None
    name = os.path.splitext(os.path.basename(filepath))[0]
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                if parts[0] == "RA":
                    ra = float(parts[1])
                elif parts[0] == "DEC":
                    dec = float(parts[1])
    if ra is None or dec is None:
        raise ValueError(f"RA/DEC not found in {filepath}")
    return [name, ra, dec]

def extract_all_from_folder(folder_path):                                # to read each file in a folder
    all_data = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".hdr"):
            filepath = os.path.join(folder_path, filename)
            data = extract_ra_dec(filepath)
            all_data.append(data)
    return all_data

#-----------DATA--------------------------------------------------------------------------------------------------
#------------CHIME------------                                           # https://www.chime-frb.ca/catalog
if catalog is not None:
    data_CHIME = catalog.as_dict()
    df_CHIME = pd.DataFrame(data_CHIME).T
else:                                                                    # in case of import problem we have saved a catalog from may 2026
    CHIME_path = os.path.join(os.path.dirname(__file__), "CHIME", "catalog.csv")
    df_CHIME = pd.read_csv(CHIME_path)

names_CHIME = df_CHIME["tns_name"].values.astype(str)                    # FRBs name
ra_CHIME = df_CHIME["ra"].values.astype(float)                           # right ascension J2000 (degrees)
raErr_CHIME = df_CHIME["ra_err"].values.astype(float)
dec_CHIME = df_CHIME["dec"].values.astype(float)                         # declination J2000 (degrees)
decErr_CHIME = df_CHIME["dec_err"].values.astype(float)
radec_CHIME = list(zip(ra_CHIME, dec_CHIME))                             # put it in a tuple

mask_rep = df_CHIME["repeater_name"] != "-9999"                          # select repeaters
radec_rep_CHIME = list(zip(ra_CHIME[mask_rep], dec_CHIME[mask_rep]))     # stock ra and dec from repeated FRBs

dm_CHIME = df_CHIME["dm_fitb"].values.astype(float)                      # DM determined using the fitting algorithm fitburst (pc/cm^3)
dm_err_CHIME = df_CHIME["dm_fitb_err"].values.astype(float)
dm_rep_CHIME = df_CHIME[mask_rep]["dm_fitb"].values.astype(float)        # for repeated FRBs
dm_nonrep_CHIME = df_CHIME[~mask_rep]["dm_fitb"].values.astype(float)    # for non-repeated FRBs (~ means the negation)

#---------PASA/FRBCAT---------                                           # https://cdsarc.cds.unistra.fr/viz-bin/cat/J/other/PASA/33.45
viz = Vizier(
    columns=["Name", "RAJ2000", "DEJ2000", "DM", "e_DM"],
    row_limit=-1
    )
table = viz.query_constraints(catalog="J/other/PASA/33.45")[0].to_pandas()

names_PASA = (table["Name"].astype(str)).to_numpy()
dm_PASA = pd.to_numeric(table["DM"], errors="coerce").to_numpy()
dm_err_PASA = pd.to_numeric(table["e_DM"], errors="coerce").to_numpy()

coords_eq = SkyCoord(                                                    # Vizier gives ra in h:M:s and dec in d:m:s >> 2 different conversion
    ra=table["RAJ2000"].values, 
    dec=table["DEJ2000"].values, 
    unit=(u.hourangle, u.deg)
    )
radec_PASA = list(zip(coords_eq.ra.deg, coords_eq.dec.deg))

warnings.filterwarnings(
    "ignore",
    message=".*OverflowError converting to IntType in column specobjid.*"
)

#------------ASKAP------------                                               # https://researchdata.edu.au/data-askap-latitude-frb-sample/1333555 
ASKAP_folder = "/Users/swan/Desktop/Cours/Theory_modules/project/data/ASKAP" # ASKAP data (FRB_XX.hdr files must be downloaded 1 by 1 and placed in this directory)

coords = extract_all_from_folder(ASKAP_folder)
names_ASKAP = [c[0] for c in coords]          
radec_ASKAP = list(zip([c[1] for c in coords], [c[2] for c in coords]))      # we stock the data needed (ie, RA-DEC in a tuple)