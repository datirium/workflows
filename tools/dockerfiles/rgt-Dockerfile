#################################################################
# Dockerfile
#
# Software:         RGT
# Software Version: v0.0.1
# Description:      Dockerized version of Regulatory Genomics Toolbox (RGT)
# Website:          https://github.com/Barski-lab/workflows
# Provides:         python   2.7.12
#                   cython   0.29.11
#                   numpy    1.16.4
#                   scipy    1.2.0
#                   pybigwig 0.3.12
#                   rgt      0.12.1
# Base Image:       biowardrobe2/scidap:v0.0.3
# Build Cmd:        docker build --rm -t biowardrobe2/rgt:v0.0.1 -f rgt-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/rgt:v0.0.1
# Run Cmd:          docker run --rm -ti biowardrobe2/rgt:v0.0.1
#################################################################


### Base Image
FROM biowardrobe2/scidap:v0.0.3
LABEL maintainer="misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive


################## BEGIN INSTALLATION ######################

WORKDIR /tmp

ENV VER_CYTHON 0.29.11
ENV VER_NUMPY 1.16.4
ENV VER_SCIPY 1.2.0
ENV VER_PYBIGWIG 0.3.12
ENV VER_RGT 0.12.1

### Installing dependencies
RUN pip install -U pip && \
    hash pip && \
    pip install setuptools && \
    pip install cython=="${VER_CYTHON}" && \
    pip install numpy=="${VER_NUMPY}" && \
    pip install scipy=="${VER_SCIPY}" && \
    pip install pyBigWig=="${VER_PYBIGWIG}" && \
### Installing RGT
    pip install rgt=="${VER_RGT}" && \
### Cleaning
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true