set -ex

cp  ../VERSION ./VERSION
cp ../environment.yaml ./environment.yaml

version=$(cat ./VERSION | sed -nre 's/^[^0-9]*(([0-9]+\.)*[0-9]+).*/\1/p')

function cleanup {
    rm  ./VERSION
    rm  ./environment.yaml
}

trap cleanup EXIT


docker build -t ceniklab/riboflow:latest . 
docker run -it ceniklab/riboflow:latest apt list | sed 's/\x1b\[[0-9;]*m//g' > ./apt.list
docker run -it ceniklab/riboflow:latest conda list > ./conda.list
docker images
