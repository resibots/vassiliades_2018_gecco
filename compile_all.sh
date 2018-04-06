#!/bin/bash
mkdir install
cd code/sferes2/exp
ln -s ../../gecco2018exp

cd ../../

INSTALL="$(realpath ../install)"
DART_PATH=$INSTALL/dart_path/
echo "Install directory: ${INSTALL}"
echo "DART path: ${DART_PATH}"

# compile dart
cd dart
mkdir build
cd build
cmake -DDART_ENABLE_SIMD=ON -DCMAKE_INSTALL_PREFIX=$DART_PATH ..
make -j4
make install

# compile hexapod_common
cd ../../hexapod_common/hexapod_controller
./waf configure --prefix=$INSTALL
./waf
./waf install

cd ../hexapod_models
./waf configure --prefix=$INSTALL
./waf
./waf install

# compile hexapod_simu
cd ../../hexapod_simu/hexapod_dart
./waf configure --prefix=$INSTALL --dart=$DART_PATH
./waf
./waf install

# compile sferes
cd ../../sferes2/
./waf configure --cpp11=yes
./waf

# compile the experiments
./waf configure --cpp11=yes --dart=$DART_PATH --exp gecco2018exp
./waf --exp gecco2018exp

cd ..
