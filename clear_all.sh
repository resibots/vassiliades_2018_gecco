#!/bin/bash
cd code/dart/build
make clean
cd ..
rm -fr build/

cd ../sferes2/
./waf clean
cd exp/
rm gecco2018exp
cd ..
git checkout -- .
git clean -df
cd ..

cd hexapod_common/hexapod_controller
./waf clean
cd ..
git checkout -- .
git clean -df
cd ..

cd hexapod_simu/hexapod_dart
./waf clean
cd ..
git checkout -- .
git clean -df
cd ..

cd ..

git checkout -- .
git clean -df
