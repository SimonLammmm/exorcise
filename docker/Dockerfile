# Exorcise: the whole toolchain in one image.
#
# Builds from the working tree, not from GitHub, so what you have checked out is
# what you get. That means the build context is the repository root and the
# Dockerfile has to be named explicitly:
#
#   docker build -f docker/Dockerfile -t simonlammmm/exorcise .
#
# For an architecture other than the host's:
#
#   docker buildx build --platform linux/amd64 -f docker/Dockerfile \
#     -t simonlammmm/exorcise:amd64 .
#
# Two conda environments, one per language. They are created in separate layers
# so that changing the R dependencies does not rebuild TensorFlow, and neither
# rebuilds when only source files change.

FROM ubuntu:24.04

SHELL ["/bin/bash", "-c"]
ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_ROOT=/opt/miniforge3
ENV EXORCISE_HOME=/exorcise

WORKDIR /

# System dependencies. curl fetches Miniforge; g++ is needed to build a couple
# of the R and Python packages from source when no wheel is available.
RUN apt-get update && \
    apt-get install --no-install-recommends -y \
      ca-certificates \
      curl \
      g++ && \
    rm -rf /var/lib/apt/lists/*

# Miniforge, which provides conda, mamba and conda-forge by default.
RUN curl -fsSL -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash "Miniforge3-$(uname)-$(uname -m).sh" -b -p "${MAMBA_ROOT}" && \
    rm -f Miniforge*.sh

# The two environments. Only the yaml files are copied at this point, so editing
# any source file leaves both of these layers cached.
COPY env/exorcise.yaml /tmp/env/exorcise.yaml
RUN "${MAMBA_ROOT}/bin/mamba" env create -f /tmp/env/exorcise.yaml -y && \
    "${MAMBA_ROOT}/bin/mamba" clean --all -y

COPY env/exorcise-py.yaml /tmp/env/exorcise-py.yaml
RUN "${MAMBA_ROOT}/bin/mamba" env create -f /tmp/env/exorcise-py.yaml -y && \
    "${MAMBA_ROOT}/bin/mamba" clean --all -y && \
    rm -rf "${MAMBA_ROOT}/pkgs" /tmp/env

# The Python environment goes on the PATH so that MAGeCK finds RRA, and the
# pipeline finds MAGeCK, without anything being patched. Earlier versions rewrote
# three source files with sed at build time to hard-code these paths.
ENV PATH="${MAMBA_ROOT}/envs/exorcise-py/bin:${PATH}"

# BLAT and twoBitToFa live in the R environment, which is not on the PATH.
# Exorcise reads these rather than being patched.
ENV EXORCISE_BLAT="${MAMBA_ROOT}/envs/exorcise/bin/blat"
ENV EXORCISE_TWOBITTOFA="${MAMBA_ROOT}/envs/exorcise/bin/twoBitToFa"

# Now the source. env/ is included so the pinned versions can be inspected from
# inside a running container, and data/ because the bundled human exome is small
# and immediately useful. example/ and manuscript/ are not copied.
COPY R/ ${EXORCISE_HOME}/R/
COPY bin/ ${EXORCISE_HOME}/bin/
COPY py/ ${EXORCISE_HOME}/py/
COPY env/ ${EXORCISE_HOME}/env/
COPY data/ ${EXORCISE_HOME}/data/
COPY docker/entrypoint.sh ${EXORCISE_HOME}/docker/entrypoint.sh
COPY README.md CHANGELOG.md LICENSE ${EXORCISE_HOME}/

RUN chmod +x "${EXORCISE_HOME}"/bin/* "${EXORCISE_HOME}/docker/entrypoint.sh" && \
    "${MAMBA_ROOT}/envs/exorcise-py/bin/pip" install --no-cache-dir "${EXORCISE_HOME}/py"

ENTRYPOINT ["/exorcise/docker/entrypoint.sh"]

# Bind mount your working directory here.
VOLUME ["/data/"]
WORKDIR /data/
