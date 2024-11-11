FROM continuumio/miniconda3

RUN conda install -c bioconda -c conda-forge \
    seqkit=2.8.2 \
    biopython=1.84 \
    fastp=0.23.4 \
    cutadapt=4.9 \
    levenshtein=0.26.0 \
    regex=2024.9.11 \
    snapgene-reader=0.1.21 \
    scipy=1.14.1 \
    pandas=2.2.2 \
    click=8.1.7





