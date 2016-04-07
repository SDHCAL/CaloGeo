# CaloGeo
# DD4hep calorimeter geometry implementations
Copyright IPNL/CNRS/IN2P3

## INSTALL:

```bash
mkdir build
cd build
cmake -C $ILCSOFT/ILCSoft.cmake ..
make install
```

#### Display geometry

Use DD4hep package binary geoDisplay :
```bash
. ./bin/thisCaloGeo.sh
geoDisplay -compact ./SDHCAL/compact/sdhcal_m3.xml
```

### Bug report

You can send emails to <rete@ipnl.in2p3.fr>
or use the [github issues interface](https://github.Com/DQM4HEP/CaloGeo/issues)
