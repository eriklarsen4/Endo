# Base image with R 4.4.1 (required for MASS >= 7.3-65)
FROM rocker/rstudio:4.4.1

# Set working directory inside the container
WORKDIR /home/rstudio/endo

# Copy the entire Endo repo into the container
COPY . /home/rstudio/endo

# Avoid interactive prompts during apt operations
ENV DEBIAN_FRONTEND=noninteractive

# Install system libraries required by ragg and other graphics packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libwebp-dev \
    libwebpmux3 \
    libpng-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install renv and restore environment
RUN Rscript -e "install.packages('renv', repos = 'https://cloud.r-project.org')" \
    && Rscript -e "renv::restore()"

# Install Python and pip
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and install the Python package from the python/ directory
RUN python3 -m pip install --upgrade pip \
    && python3 -m pip install -e python \
    && python3 -m pip install pytest pytest-cov seaborn

# Install the R package from the Rpkg directory
RUN Rscript -e "remotes::install_local('Rpkg', upgrade = 'never')"

# Expose RStudio Server port
EXPOSE 8787

# Ensure entrypoint script is executable
RUN chmod +x /home/rstudio/endo/entrypoint.sh

# Use the entrypoint script
ENTRYPOINT ["/home/rstudio/endo/entrypoint.sh"]