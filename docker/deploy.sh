
docker login -u hakanozadam

version=$(cat ../VERSION | sed -nre 's/^[^0-9]*(([0-9]+\.)*[0-9]+).*/\1/p')
echo "version: $version"

# push the image
docker push hakanozadam/riboflow:latest
docker push hakanozadam/riboflow:$version
