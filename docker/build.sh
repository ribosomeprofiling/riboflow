#!/bin/bash
set -ex

# Run from the repo root regardless of where the script is invoked from
cd "$(dirname "$0")/.."

# Extract version
version=$(cat ./VERSION | sed -nre 's/^[^0-9]*(([0-9]+\.)*[0-9]+).*/\1/p')

# Build the image (Forcing linux/amd64 for compatibility). Tag it as BOTH
# :latest and :<VERSION>: the docker/apptainer profiles in nextflow.config pull
# :<VERSION>, so a local-only build that is not pushed via deploy.sh would
# otherwise leave the pipeline running the previously pulled image.
docker build --platform linux/amd64 \
    -t danielnguyener/riboflow:latest \
    -t danielnguyener/riboflow:${version} \
    -f docker/Dockerfile .

# Generate lists (optional, but useful)
docker run --rm danielnguyener/riboflow:latest apt list --installed 2>/dev/null > ./docker/apt.list
docker run --rm danielnguyener/riboflow:latest bash -c "source activate riboflow && conda list" > ./docker/conda.list

echo "Build complete. Tagged danielnguyener/riboflow:latest and danielnguyener/riboflow:${version}"
