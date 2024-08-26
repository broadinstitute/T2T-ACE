# Use base image with miniconda3 installed
FROM continuumio/miniconda3
MAINTAINER Yueyao Gao <tag@broadinstitute.org>

# Specify Workdir
WORKDIR /BaseImage

# Create the environment and copy the scripts
RUN mkdir -p /BaseImage/T2T-ACE/
COPY T2T_ACE/* /BaseImage/T2T-ACE/
RUN conda env create -f /BaseImage/T2T-ACE/T2T_ACE_env_explicit.yml

# Set the working directory to the T2T-ACE directory
WORKDIR /BaseImage/T2T-ACE