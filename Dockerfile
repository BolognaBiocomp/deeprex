# Base Image
FROM continuumio/miniconda3

# Metadata
LABEL base.image="continuumio/miniconda3"
LABEL version="1.0"
LABEL software="DeepREx"
LABEL software.version="2018012"
LABEL description="an open source software tool to predict residue solvent exposure from protein sequence"
LABEL website="https://deeprex.biocomp.unibo.it"
LABEL documentation="https://deeprex.biocomp.unibo.it"
LABEL license="GNU GENERAL PUBLIC LICENSE Version 3"
LABEL tags="Proteomics"
LABEL maintainer="Castrense Savojardo <castrense.savojardo2@unibo.it>"

ENV PYTHONDONTWRITEBYTECODE=true

WORKDIR /usr/src/deeprex

COPY . .

WORKDIR /data/

RUN conda update -n base conda && \
   conda install --yes nomkl keras biopython && \
   conda install -c conda-forge -c bioconda hhsuite && \
   conda clean -afy \
   && find /opt/conda/ -follow -type f -name '*.a' -delete \
   && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
   && find /opt/conda/ -follow -type f -name '*.js.map' -delete

# Verbosity level of Tensorflow
ENV TF_CPP_MIN_LOG_LEVEL=3 DEEPREX_ROOT=/usr/src/deeprex PATH=/usr/src/deeprex:$PATH

ENTRYPOINT ["/usr/src/deeprex/deeprex.py"]
