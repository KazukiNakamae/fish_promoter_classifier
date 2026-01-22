FROM continuumio/miniconda3:25.3.1-1

RUN printf "channels:\n  - conda-forge\n  - bioconda\nchannel_priority: strict\n" > /opt/conda/.condarc

RUN conda create -y -n dnabert-env \
      python=3.10 \
      meme=5.5.9 \
      bedtools=2.31.1 \
      parallel \
      seqkit=2.12.0 \
    && conda clean --all --yes

SHELL ["bash", "-lc"]
ENV PATH=/opt/conda/envs/dnabert-env/bin:$PATH
