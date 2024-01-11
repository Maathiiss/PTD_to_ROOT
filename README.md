# Falaise Skeleton Modules

Example of Falaise modules retrieving and browsing different data banks: SD, CD, etc ...

`falaise-skeleton-module-sd.cc`
Skeleton module browsing SD bank (Simulated Data)

`falaise-skeleton-module-cd.cc`
Skeleton module browsing CD bank (Calibrated Data)

`falaise-skeleton-module-tcd.cc`
Skeleton module browsing TCD bank (Tracker Clustering Data)

`falaise-skeleton-module-ttd.cc`
Skeleton module browsing TTD bank (Tracker Trajectory Data)

`falaise-skeleton-module-ptd.cc`
Skeleton module browsing PTD bank (Particle Track Data)

Other modules to come soon for UDD and pCD banks ...


Compilation and usage at CC/IN2P3 using Falaise 5.1:
```
source setup_ccin2p3_falaise51.sh

mkdir build-falaise51
cd build-falaise51/
cmake ../
make
cd ../

# run the demo pipeline with a nice 2b0n simulated event
flreconstruct -p build-falaise51/falaise-skeleton-module-pipeline.conf -i brio/mc-se82-2b0n-event.brio
```
