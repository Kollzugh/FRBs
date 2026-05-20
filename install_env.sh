# !!See README.md to install the setup with this file!!

#!/usr/bin/env bash
set -euo pipefail

# Initialization of the Conda environment
read -rp "Name of the Conda environment to create: " ENV_NAME
if [ -z "$ENV_NAME" ]; then
    echo "Error: the environment name cannot be empty."
    exit 1
fi
conda create --name "$ENV_NAME" python=3.9.23 -y
source "$(conda info --base)/etc/profile.d/conda.sh"

# Install environment-specific dependencies
conda activate "$ENV_NAME"
python -m pip install --no-build-isolation git+https://github.com/desihub/desiutil.git@3.2.5
python -m pip install --no-build-isolation -r requirements.txt
