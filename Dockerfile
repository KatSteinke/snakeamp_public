FROM condaforge/mambaforge:latest

# TODO: change to use separate container

RUN mkdir /conda-envs
COPY envs/ /conda-envs

RUN mamba install -c defaults -c bioconda -c conda-forge -c anaconda -n base pandas<3.* pyarrow openpyxl snakemake==7.22.0 numpy>=1.22.2 && \
    mamba env create --file /conda-envs/emu_env.yml && \
    mamba env create --file /conda-envs/nanopore_qc.yml && \
	mamba clean --all -y

WORKDIR /snake_data
COPY . /snake_data/