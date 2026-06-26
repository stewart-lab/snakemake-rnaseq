## base image with Debian + Python 3.12
FROM condaforge/miniforge3:latest

## set strict channel priorities for conda
RUN conda config --set channel_priority strict

## install snakemake
RUN mamba install -y -c conda-forge -c bioconda snakemake graphviz

# install chrome dependencies (for kaleido/plotly to export plot .png/.svg)
RUN apt update && apt-get install -y \
    libnss3 \
    libatk-bridge2.0-0 \
    libcups2 \
    libxcomposite1 \
    libxdamage1 \
    libxfixes3 \
    libxrandr2 \
    libgbm1 \
    libxkbcommon0 \
    libpango-1.0-0 \
    libcairo2

# install python packages
RUN pip install --upgrade pip
RUN pip install kaleido

# install chrome for kaleido
RUN kaleido_get_chrome

## copy the pipeline code into the container
WORKDIR /app
COPY app /app

ENTRYPOINT ["python", "main.py"]