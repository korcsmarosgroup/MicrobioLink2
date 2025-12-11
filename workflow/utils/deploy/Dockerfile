
# Dockerfile for MicrobioLink2 Preprocessing Module
# Installs Python and R dependencies for all scripts in workflow/utils/preprocessing

FROM rocker/r-base:4.3.2

# System dependencies for Python, R, and common bioinformatics tools
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 python3-pip python3-venv python3-dev \
        build-essential libxml2-dev libssl-dev libcurl4-openssl-dev \
        git wget unzip locales \
        libhdf5-dev \
        libz-dev \
        libpng-dev \
        libblas-dev liblapack-dev gfortran \
        r-cran-yaml r-cran-matrix r-cran-readr r-cran-processx \
        vim-common \
    && rm -rf /var/lib/apt/lists/*

# Set UTF-8 locale
RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
	locale-gen
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8


# Create and activate a Python virtual environment
RUN python3 -m venv /venv
ENV PATH="/venv/bin:$PATH"

# Install Python dependencies in the venv
RUN pip install --no-cache-dir \
    scanpy \
    matplotlib \
    pandas \
    numpy \
    scipy \
    pyyaml \
    scrublet

# Alternatively, get STAR source using git
RUN git clone https://github.com/alexdobin/STAR.git && \
    cd STAR/source && \
    make STAR
ENV PATH="/STAR/source:$PATH"

# Install R dependencies using renv if lockfile exists
RUN if [ -f Rpreprocessing/renv.lock ]; then \
	  R -e 'install.packages("renv", repos="https://cloud.r-project.org/"); renv::restore(lockfile = "Rpreprocessing/renv.lock")'; \
	fi

# Install additional R packages from DESCRIPTION (if not in renv)
RUN R -e 'install.packages(c("Matrix", "processx", "readr", "Seurat", "this.path", "tools", "yaml"), repos="https://cloud.r-project.org/")'

# Copy preprocessing scripts
COPY preprocessing /preprocessing
WORKDIR /preprocessing
