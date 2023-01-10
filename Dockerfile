FROM condaforge/mambaforge:latest

RUN mamba install -c defaults -c bioconda -c conda-forge -n base emu=3.4.4 filtlong=0.2.1 seqtk=1.3

COPY emu_barcodes.sh emu_barcodes.sh