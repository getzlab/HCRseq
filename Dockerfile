# 1. Use a lightweight Mamba base image for faster, more reliable builds
FROM mambaorg/micromamba:1.5.8-bookworm-slim

USER root

# 2. Set the working directory
WORKDIR /usr/src/app

# 3. Copy the environment file and create the Conda environment.
#    This single step replaces all manual compilation and the pip requirements file.
#COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yaml .
COPY environment.yaml .
RUN micromamba install -y -n base -f environment.yaml && \
    micromamba clean --all --yes

# 4. Activate the Conda environment for all subsequent commands
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# 5. Copy your application code into the container
#COPY --chown=$MAMBA_USER:$MAMBA_USER hcrseq ./hcrseq
#COPY --chown=$MAMBA_USER:$MAMBA_USER setup.py .
COPY hcrseq ./hcrseq
COPY setup.py .

# 6. Install your local package into the Conda environment.
#    --no-deps is important because Conda already handled all dependencies.
RUN pip install --no-deps -e .

ENV PATH=/opt/conda/bin:$PATH

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh" ]
CMD [ "/bin/bash" ]
