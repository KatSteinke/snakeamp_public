FROM condaforge/miniforge3:26.9.4-0

# TODO: change to use separate container

RUN mkdir /conda-envs
COPY envs/ /conda-envs

RUN mamba install -c bioconda -c conda-forge -n base "ca-certificates>=2026.10.19" "certifi>=2026.10.19" pandas=2.* "pip>=26.2.2" pyarrow openpyxl snakemake==9.27.1 "numpy>=1.22.2" && \
    mamba env create --file /conda-envs/emu_env.yml && \
    mamba env create --file /conda-envs/kraken_env.yml && \
    mamba env create --file /conda-envs/nanopore_qc.yml && \
	mamba clean --all -y

WORKDIR /snake_data
COPY . /snake_data/