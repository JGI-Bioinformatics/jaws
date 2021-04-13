FROM continuumio/miniconda3:4.8.2
MAINTAINER po-e@lanl.gov

LABEL version="1.0.0"
LABEL software="nmdc_taxa_profilers"
LABEL tags="bioinformatics"

ENV container docker

RUN apt-get update -y \
    && apt-get install -y build-essential unzip wget curl gawk \
    && apt-get clean

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

# install gottcha2
RUN conda install minimap2 pandas
RUN wget https://github.com/poeli/GOTTCHA2/archive/2.1.7.tar.gz \
    && tar -xzf 2.1.7.tar.gz \
    && cp GOTTCHA2-2.1.7/*.py /usr/local/bin \
    && rm -rf GOTTCHA2-2.1.7/ 2.1.7.zip

# install kraken2
RUN conda install kraken2=2.1.0

# install centrifuge
RUN wget https://github.com/DaehwanKimLab/centrifuge/archive/v1.0.4-beta.tar.gz \
    && tar -xzf v1.0.4-beta.tar.gz \
    && cd centrifuge-1.0.4-beta \
    && make install prefix=/usr/local

# install krona
RUN conda install krona \
    && ktUpdateTaxonomy.sh

CMD ["/bin/bash"]
