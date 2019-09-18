#!bin/sh

cd src
make clean
make all
../bin/dEdx.exe -b -i /eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/2019_06_14/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_v2.root -s -r
../bin/dEdx.exe -b -i $DESY_out/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_v2.root -r

../bin/SpatialResol.exe -b -t0 -i /eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/2019_06_14/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_v2.root -s -r
../bin/SpatialResol.exe -b -t0 -i $DESY_out/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_v2.root -r

../bin/SpatialResol.exe -b -t0 -i /eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/2019_06_14/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_v2.root -r

../bin/SpatialResol.exe -b -t0 -i /eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/2019_06_14/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/v6/debug_iter0.root -r
../bin/SpatialResol.exe -b -t1 -i /eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/2019_06_14/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_v2_iter1.root -r
../bin/SpatialResol.exe -b -t2 -i /eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/2019_06_14/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_v2_iter2.root -r
../bin/SpatialResol.exe -b -t3 -i /eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/2019_06_14/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_v2_iter3.root -r
../bin/SpatialResol.exe -b -t4 -i /eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/2019_06_14/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_v2_iter4.root -r

../bin/SpatialResol.exe -b -t0 -i $DESY_out/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_iter0.root -r
../bin/SpatialResol.exe -b -t1 -i $DESY_out/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_iter1.root -r
../bin/SpatialResol.exe -b -t2 -i $DESY_out/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_iter2.root -r
../bin/SpatialResol.exe -b -t3 -i $DESY_out/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_iter3.root -r
../bin/SpatialResol.exe -b -t4 -i $DESY_out/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/debug_iter4.root -r

../bin/SpatialResol.exe -b -t0 -i $DESY_out/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/test_iter0.root -r -c
../bin/SpatialResol.exe -b -t1 -i $DESY_out/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/test_iter1.root -r -c
../bin/SpatialResol.exe -b -t2 -i $DESY_out/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/test_iter2.root -r -c
../bin/SpatialResol.exe -b -t3 -i $DESY_out/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/test_iter3.root -r -c
../bin/SpatialResol.exe -b -t4 -i $DESY_out/R2019_06_14-16_42_26-000.root -o /eos/user/s/ssuvorov/DESY_testbeam/test_iter4.root -r -c


