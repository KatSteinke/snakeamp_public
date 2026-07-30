FROM condaforge/miniforge3:latest

# TODO: change to use separate container

RUN mkdir /conda-envs
COPY envs/ /conda-envs

RUN mamba install -c bioconda -c conda-forge -n base pandas=2.* pyarrow openpyxl snakemake==9.14 "numpy>=1.22.2" && \
    mamba env create --file /conda-envs/emu_env.yml && \
    mamba env create --file /conda-envs/kraken_env.yml && \
    mamba env create --file /conda-envs/nanopore_qc.yml && \
	mamba clean --all -y

WORKDIR /snake_data
COPY . /snake_data/