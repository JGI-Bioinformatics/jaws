FROM continuumio/miniconda2

# install software
RUN conda install -c bioconda bbmap
RUN conda install -c bioconda samtools

# move data and script into container
WORKDIR /bbmap
COPY script.sh /usr/local/bin/script.sh


