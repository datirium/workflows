#################################################################
# Dockerfile
#
# Software:         Python Pandas
# Software Version: v0.0.1
# Description:      Dockerized version of Python Pandas Package
# Website:          https://github.com/Barski-lab/workflows
# Provides:         Python 2.7.12
#                   Pandas==0.23.4
# Base Image:       biowardrobe2/scidap:v0.0.3
# Build Cmd:        docker build --rm -t biowardrobe2/python-pandas:v0.0.1 -f python-pandas-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/python-pandas:v0.0.1
# Run Cmd:          docker run --rm -ti biowardrobe2/python-pandas:v0.0.1
#################################################################


### Base Image
FROM biowardrobe2/scidap:v0.0.3
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive


################## BEGIN INSTALLATION ######################

WORKDIR /tmp

ENV VER_PANDAS 0.23.4

### Installing Pandas
RUN pip install pandas=="${VER_PANDAS}" && \

### Cleaning
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true
