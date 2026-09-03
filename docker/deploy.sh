#!/bin/bash

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

IMAGE_NAME="riboflow"

while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--name)
            IMAGE_NAME="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  -n, --name NAME    Image name (default: riboflow)"
            echo "                     Example: myusername/riboflow"
            echo "  -h, --help         Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

VERSION=$(grep "version =" "$PROJECT_ROOT/VERSION" | cut -d"'" -f2)

if [ -z "$VERSION" ]; then
    echo "Error: Could not extract version from $PROJECT_ROOT/VERSION"
    exit 1
fi

echo "Deploying Image: $IMAGE_NAME"
echo "Version: $VERSION"
echo ""

if [[ "$(docker images -q danielnguyener/riboflow:latest 2> /dev/null)" == "" ]]; then
    echo "Error: Local image 'danielnguyener/riboflow:latest' not found."
    echo "Please run ./docker/build.sh first."
    exit 1
fi

echo "Tagging local image 'danielnguyener/riboflow:latest' as '${IMAGE_NAME}:latest'..."
docker tag danielnguyener/riboflow:latest "${IMAGE_NAME}:latest"

echo "Tagging local image 'danielnguyener/riboflow:latest' as '${IMAGE_NAME}:${VERSION}'..."
docker tag danielnguyener/riboflow:latest "${IMAGE_NAME}:${VERSION}"

echo ""
echo "Pushing 'latest' tag..."
docker push "${IMAGE_NAME}:latest"

echo "Pushing '${VERSION}' tag..."
docker push "${IMAGE_NAME}:${VERSION}"

echo ""
echo "Deployment complete!"
