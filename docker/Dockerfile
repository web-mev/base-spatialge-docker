FROM mambaorg/micromamba

# for building the conda env
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
# Build the environment based on the yml file
RUN micromamba install -y -n base -f /tmp/environment.yml \
  && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN R -e 'install.packages("spaMM", repos="http://cran.us.r-project.org")'
RUN R -e 'library("devtools"); devtools::install_github("fridleylab/spatialGE", quiet=TRUE, dependencies=TRUE)'

ADD prep_stlist.R /usr/local/bin/
