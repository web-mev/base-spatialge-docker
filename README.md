# base-spatialge-docker

This repository contains a build for a Docker image that contains all installations and dependencies for the spatialGE tool (https://github.com/FridleyLab/spatialGE).

Note that the micromamba container from which this image derives has an implicit entrypoint command which activates the conda environment contained in the image (see https://micromamba-docker.readthedocs.io/en/latest/quick_start.html#activating-a-conda-environment-for-entrypoint-commands). If you invoke the image and override this entrypoint, the conda environment will not be activated and used.