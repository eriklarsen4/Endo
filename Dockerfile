# Base image with RStudio Server and R 4.4.1
FROM rocker/rstudio:4.4.1

# Avoid interactive prompts during apt operations
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies for tidyverse, renv, and Python
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-venv \
        git \
        make \
        g++ \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libxt-dev \
        libpng-dev \
        libjpeg-dev \
        libtiff5-dev \
        libfontconfig1-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create Python virtual environment for the Python side of the Endo package
ENV VIRTUAL_ENV=/home/rstudio/pyenv

RUN python3 -m venv "${VIRTUAL_ENV}" && \
    . "${VIRTUAL_ENV}/bin/activate" && \
    pip install --upgrade pip

# Ensure the venv is on PATH
ENV PATH="${VIRTUAL_ENV}/bin:${PATH}"

# Set working directory inside the container
WORKDIR /home/rstudio/endo

# Copy the entire Endo repo into the container
COPY . /home/rstudio/endo

# System libraries required by ragg and other graphics packages
RUN apt-get update && apt-get install -y \
    libwebp-dev \
    libwebpmux3 \
    libpng-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libfontconfig1-dev

# Install renv and restore environment
RUN Rscript -e "install.packages('renv', repos = 'https://cloud.r-project.org')" && \
    Rscript -e "renv::restore()"


# Install the local R package (the repo root is the package)
RUN Rscript -e "install.packages('remotes', repos = 'https://cloud.r-project.org')" && \
    Rscript -e "remotes::install_local('.', upgrade = 'never')"

# Install the Python package if present
RUN if [ -f "pyproject.toml" ] || [ -f "setup.cfg" ] || [ -f "setup.py" ]; then \
        pip install -e . ; \
    fi

# Copy entrypoint script into the container
COPY entrypoint.sh /home/rstudio/endo/entrypoint.sh
RUN chmod +x /home/rstudio/endo/entrypoint.sh

# Set entrypoint to activate environments and run startup checks
ENTRYPOINT ["/home/rstudio/endo/entrypoint.sh"]


# Ensure the rstudio user owns the project directory
RUN chown -R rstudio:rstudio /home/rstudio/endo

# Healthcheck to verify RStudio Server is responding
HEALTHCHECK CMD curl --fail http://localhost:8787 || exit 1

# Switch to the rstudio user for interactive work
USER rstudio
WORKDIR /home/rstudio/endo

# Default command opens a shell
CMD ["/bin/bash"]