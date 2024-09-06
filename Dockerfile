# Use base image with miniconda3 installed
FROM continuumio/miniconda3
MAINTAINER Yueyao Gao <tag@broadinstitute.org>

# Specify Workdir
WORKDIR /BaseImage

# Create the environment and copy the scripts
RUN mkdir -p /BaseImage/T2T-ACE/ACElib/
COPY T2T_ACE/ACElib/ /BaseImage/T2T-ACE/ACElib/
COPY T2T_ACE/T2T_ACE_env_explicit.yml /BaseImage/T2T-ACE/T2T_ACE_env_explicit.yml
COPY T2T_ACE/run_T2T-ACE.py /BaseImage/T2T-ACE/run_T2T-ACE.py
RUN conda env create -f /BaseImage/T2T-ACE/T2T_ACE_env_explicit.yml

# Set the working directory to the T2T-ACE directory
WORKDIR /BaseImage/T2T-ACE