# Falaise Skeleton Modules

Example of Falaise modules retrieving and browsing data banks: SD, CD, etc ...

Compilation and usage at CC/IN2P3 using Falaise 5.1:
```
source setup_ccin2p3_falaise51.sh

mkdir build-falaise51
cd build-falaise51
cmake ..
make
cd ..

flreconstruct -p build-falaise51/falaise-skeleton-module-pipeline.conf -i INPUT_FILE.brio
```
