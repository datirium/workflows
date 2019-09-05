#################################################################
# Dockerfile
#
# Software:         CLIPper
# Software Version: v0.0.1
# Description:      Dockerized version of CLIPper Python Package
# Website:          https://github.com/yeolab/clipper
# Provides:         CLIPper (latest commit 997fe25532a5bdcf5957f2a467ca283bbd550303)
# Base Image:       biowardrobe2/clipper:v0.0.1
# Build Cmd:        docker build --rm -t biowardrobe2/clipper:v0.0.1 -f clipper-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/clipper:v0.0.1
# Run Cmd:          docker run --rm -ti biowardrobe2/clipper:v0.0.1
#################################################################


### Base Image
FROM biowardrobe2/scidap:v0.0.3
LABEL maintainer="misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive


ENV VER_CLIPPER 997fe25532a5bdcf5957f2a467ca283bbd550303
ENV URL_CLIPPER git+git://github.com/YeoLab/clipper.git@${VER_CLIPPER}


################## BEGIN INSTALLATION ######################
WORKDIR /tmp

### Installing dependencies
RUN apt-get clean all &&\
    apt-get update &&\
    apt-get install -y libbz2-dev liblzma-dev python-tk && \
    pip install -U pip && \
    hash pip && \
    pip install setuptools numpy cython pandas && \
### Installing CLIPper
    pip install ${URL_CLIPPER} && \
### Update scikit-learn to proper (old) version
    pip install -U scikit-learn==0.17.1 && \
### Cleaning
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true
