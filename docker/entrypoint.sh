#!/bin/bash
# Activate Python virtual environment
source /home/rstudio/pyenv/bin/activate

# Optional: verify R package loads
Rscript -e "library(Endo)"

# Execute whatever command was passed (RStudio Server starts automatically)
exec "$@"

