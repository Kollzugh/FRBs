# Cosmology with Fast Radio Bursts
Author: Swan Dubief, MSc Theoretical Physics, University of Geneva, 2026.  
Context: Theory modules.  

## Installation/Setup on Unix/Linux systems with conda:
### Installation  
```
chmod +x install_env.sh
./install_env.sh
````

### Initialization  
Activate the environment to use codes of this project:
```
conda activate <env_name>
```  
Finally, go into the following directory: /project

## General content
In the project directory there are:
- this README.file (to launch it on VSCode: cmd+shift+V).
- a gitignore file to ensure that certain files not tracked by Git remain untracked.
- a data directory that contains data from ASKAP, a CHIME catalog to ensure we have data in case of authorization/server problems, an exemple of the original table from CHIME, the data.py file to import FRB data, and the FRBs_galaxies_catalog.txt file which is the cross-match catalog used for the best fit model. Note you can generate your own cross-match catalog (see examples below).
- an outputs directory that contains useful example graphs, the beamer presentation and the report.
- the install_env.sh and the requirements.txt file both used to install the basic setup (see installation above).
- 5 .py files (utils.py get a lot of useful functions, constants.py contains all the constants from this project, mockFRB.py to create a simulated catalog of FRBs, plotDM.py to plot the theoretical Macquart relation in different ways, plotREALDATA.py to plot different graphs of real data (Mollweide, histogram, DM(z)) and also to create a cross-match catalog of FRBs/hosts).


## Exemples  
- mockFRB.py:  
python mockFRB.py  

- plotDM.py:  
python plotDM.py --model MP18 --type cosmic  
python plotDM.py --model MP18 --type obs  
python plotDM.py --model MP18_var  

- plotREALDATA.py:  
python plotREALDATA.py --plot mollweide --coord galactic  
python plotREALDATA.py --plot mollweide --coord equatorial  
python plotREALDATA.py --plot DM_histogram  
python plotREALDATA.py --plot DM_z  (you can plot variance-z and Omegab/fIGM-H0 but you have to uncoment specific parts in red in this file)  
python plotREALDATA.py --plot DM_z --catalog yes (note you can change the tunable parameters in constants.py for the cross-match)  

!You will probably have a warning: Unable to download CHIME/FRB Catalog. We solve this by adding a local CHIME catalog (version may 2026).!