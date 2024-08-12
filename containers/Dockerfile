# Use the miniconda container
FROM condaforge/mambaforge:latest as cracofunew

LABEL description = "Container for running CELEBI python scripts"

# Install curl
RUN apt-get -y update

# set environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# copy project requirements
COPY ./cracofunew.yml .

# Create environment and activate it
RUN mamba env update --name base --file cracofunew.yml