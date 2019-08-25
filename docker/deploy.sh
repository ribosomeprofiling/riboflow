
docker login -u ceniklab

version=$(cat ../VERSION | sed -nre 's/^[^0-9]*(([0-9]+\.)*[0-9]+).*/\1/p')
echo "version: $version"

# push the image
docker push ceniklab/riboflow:latest
docker push ceniklab/riboflow:$version
