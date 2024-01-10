# Base image
# Can probably be migrated to a smaller image
FROM ubuntu:jammy-20231211.1

# Install all software under /opt/software:
RUN mkdir -p /opt/software 

# Get the miniconda installer script and run it in silent mode:
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/software/miniconda.sh && \
  bash /opt/software/miniconda.sh -b -p /opt/software/miniconda
ENV PATH="/opt/software/miniconda/bin:${PATH}"

# Install DESeq2 and dependencies:
ADD build.R /opt/software/
ADD environment.yml /opt/software/

# Build the deseq2 environment based on the yml file
RUN conda env create -f /opt/software/environment.yml

# Add to the PATH. Otherwise, conda can swallow return codes
ENV PATH="/opt/software/miniconda/envs/deseq2/bin:${PATH}"

# Open conda in shell to activate env to install spaMM, 
# then devtools-build spatialGE
SHELL ["conda", "run", "-n", "spatialge", "/bin/bash", "-c"]
RUN Rscript build.R