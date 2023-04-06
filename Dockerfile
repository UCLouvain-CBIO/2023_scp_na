## BIOCONDUCTOR
## Pull the development version of the Bioconductor Docker image
FROM bioconductor/bioconductor_docker:RELEASE_3_17

RUN apt-get update \
	## Remove packages in '/var/cache/' and 'var/lib'
	## to remove side-effects of apt-get update
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

## Change home directory
WORKDIR /home/rstudio/

COPY install_dependencies.R install_dependencies.R
COPY utils.R utils.R

## Install R dependencies
RUN Rscript install_dependencies.R
## Do a second pass in case some packages failed to installed
RUN Rscript install_dependencies.R
