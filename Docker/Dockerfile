FROM mambaorg/micromamba:2.5-debian13-slim

ARG TARGETARCH
ENV TARGETARCH=${TARGETARCH}
ENV MAMBA_ROOT_PREFIX=/home/mamba
ENV ENV_NAME=bip

#set working directory
WORKDIR /BIP

# Copy files from this directory into image
COPY . .

#mamba environment
RUN micromamba create -y -n bip -c conda-forge -c bioconda -c anaconda \
    python=3.11 \
    snakemake-minimal=9.16.3 \
    r-base=4.5.2 \
    pear=0.9.6 \
    seqkit=2.12.0 \
    libgcc=15.2.0 \
    cutadapt=5.2 \
    seqtk=1.5 \
    vsearch=2.30.4 \
    parallel=20260122 \
    rename=1.601 \
    bc=1.07.1 \
    zip=3.0 \
    chopper=0.13.0 \
# and python libs.
    pandas=3.0.0 \
    openpyxl=3.1.5 \
    biopython=1.86 \
# and all the Rpkgs
    r-tidyr=1.3.2 \
    r-data.table=1.17.8 \
    bioconductor-biostrings=2.78 \ 
    r-scales=1.4.0 \
    r-ggplot2=4.0.2 \
    r-patchwork=1.3.2 \
    r-openxlsx=4.2.8.1 \
    r-plyr=1.8.9 \
    r-readxl=1.4.5 \ 
    r-plotly=4.12.0 \
    r-forcats=1.0.1 \
    r-shiny=1.13.0 \
    r-kableextra=1.4.0 \
    r-treemapify=2.6.0 \
    r-paletteer=1.7.0 \
    r-htmltools=0.5.9 \
    r-dt=0.34.0 \
    r-leaflet=2.2.3 \
    conda-forge::c-compiler \
    conda-forge::cxx-compiler \
    bioconductor-decipher=3.6.0 \
    bioconductor-msa=1.42.0 \
    conda-forge::r-readr=2.2.0 \
    r-cowplot=1.2.0 \
    r-gridextra=2.3.1 \
    r-gtable=0.3.6 \
    && micromamba clean --all --yes

# bioconductor-muscle has a bioconda build for linux-64 (amd64) but not for
# linux-aarch64 (arm64). On amd64 install the prebuilt binary directly; on
# arm64 compile it from source using the c-compiler/cxx-compiler toolchain
# pulled in above.
RUN if [ "$TARGETARCH" = "amd64" ]; then \
        echo "Installing bioconductor-muscle from bioconda (amd64)"; \
        micromamba install -y -n bip -c conda-forge -c bioconda -c anaconda bioconductor-muscle && \
        micromamba clean --all --yes; \
    else \
        echo "Compiling muscle from source ($TARGETARCH)"; \
        micromamba run -n bip Rscript -e ' \
            if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org"); \
            BiocManager::install("muscle", type = "source", update = FALSE, ask = FALSE)'; \
    fi

USER root


#Activate BIP environment.
RUN echo "micromamba activate bip" >> ~/.bashrc

#Make directory structure
RUN mkdir SCRIPTS
RUN mkdir Barcoding

#Move files to correct place
RUN mv *.sh SCRIPTS
RUN mv *.py SCRIPTS
RUN mv *.R SCRIPTS
RUN mv *.fastq.gz Barcoding
RUN mv *.xlsx Barcoding

# Ensure BIP environment when image is run.
SHELL ["micromamba", "run", "-n", "bip", "/bin/bash", "-c"]

#Download and unzip latest BOLDDistilled lib (curl needs map environment running, hence placement here at the end)
RUN curl https://us-sea-1.linodeobjects.com/boldistilled/sintax.zip --output reflib_latest.zip
RUN python -m zipfile -e reflib_latest.zip .

#change dir name to REFS, remove zip file.
RUN mv sintax/ REFS/
RUN python -m zipfile -e reference_seqs_327K.fasta.zip ./REFS/
RUN rm reflib_latest.zip reference_seqs_327K.fasta.zip