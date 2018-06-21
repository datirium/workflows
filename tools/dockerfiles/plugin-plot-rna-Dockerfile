#################################################################
# Dockerfile
#
# Software:         Plot RNA R scripts
# Software Version: v0.0.2
# Description:      R script image for BioWardrobe plugins to plot RNA-Seq results
# Website:          https://github.com/Barski-lab/workflows
# Provides:         plugin_plot_rna.R script
# Base Image:       biowardrobe2/scidap-addons:v0.0.6
# Build Cmd:        docker build --rm -t biowardrobe2/plugin-plot-rna:v0.0.2 -f plugin-plot-rna-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/plugin-plot-rna:v0.0.2
# Run Cmd:          docker run --rm -ti biowardrobe2/plugin-plot-rna:v0.0.2 plot_rna.R --help
#################################################################


#################################################################
# Updates
#
# v0.0.1 - Generate plots, export data to text files
# v0.0.2 - Export correct RPKM distribution histogram data
#################################################################


### Base Image
FROM biowardrobe2/scidap-addons:v0.0.6
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive


################## BEGIN INSTALLATION ######################

WORKDIR /tmp

COPY ./scripts/plugin_plot_rna.R /usr/local/bin/plot_rna.R

### Installing plot_rna.R R script
RUN chmod +x /usr/local/bin/plot_rna.R && \

### Cleaning
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true
