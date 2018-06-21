#################################################################
# Dockerfile
#
# Software:         Plot DNA R scripts
# Software Version: v0.0.1
# Description:      R script image for BioWardrobe plugins to plot ChIP-Seq results
# Website:          https://github.com/Barski-lab/workflows
# Provides:         plugin_plot_dna.R script
# Base Image:       biowardrobe2/scidap-addons:v0.0.6
# Build Cmd:        docker build --rm -t biowardrobe2/plugin-plot-dna:v0.0.1 -f plugin-plot-dna-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/plugin-plot-dna:v0.0.1
# Run Cmd:          docker run --rm -ti biowardrobe2/plugin-plot-dna:v0.0.1 plot_dna.R --help
#################################################################


### Base Image
FROM biowardrobe2/scidap-addons:v0.0.6
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive


################## BEGIN INSTALLATION ######################

WORKDIR /tmp

COPY ./scripts/plugin_plot_dna.R /usr/local/bin/plot_dna.R

### Installing plot_dna.R R script
RUN chmod +x /usr/local/bin/plot_dna.R && \

### Cleaning
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true
