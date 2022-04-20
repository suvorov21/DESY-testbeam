set ver v10

set sample z10
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter0.root -r -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter1.root -r -b

set sample z20
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter0.root -r -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter1.root -r -b

set sample z30
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter0.root -r -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter1.root -r -b

set sample z40
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter0.root -r -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter1.root -r -b

set sample z50
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter0.root -r -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter1.root -r -b

set sample z60
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter0.root -r -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter1.root -r -b

set sample z70
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter0.root -r -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter1.root -r -b

set sample z80
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter0.root -r -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter1.root -r -b

set sample z90
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter0.root -r -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter1.root -r -b

set sample z95
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter0.root -r -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_{$sample}_{$ver}.root -o ~/DATA/{$sample}_{$ver}_iter1.root -r -b


./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_90deg_{$ver}.root -o ~/DATA/90deg_{$ver}_iter0.root -r --param ../params/invert_slope.ini -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_90deg_{$ver}.root -o ~/DATA/90deg_{$ver}_iter1.root -r --param ../params/invert_slope.ini -b

./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_80deg_{$ver}.root -o ~/DATA/80deg_{$ver}_iter0.root -r --param ../params/invert_slope.ini -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_80deg_{$ver}.root -o ~/DATA/80deg_{$ver}_iter1.root -r --param ../params/invert_slope.ini -b

./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_70deg_{$ver}.root -o ~/DATA/70deg_c_{$ver}_iter0.root -r --param ../params/3by1.ini -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_70deg_{$ver}.root -o ~/DATA/70deg_c_{$ver}_iter1.root -r --param ../params/3by1.ini -b
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_70deg_{$ver}.root -o ~/DATA/70deg_d_{$ver}_iter0.root -r --param ../params/diag.ini -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_70deg_{$ver}.root -o ~/DATA/70deg_d_{$ver}_iter1.root -r --param ../params/diag.ini -b
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_70deg_{$ver}.root -o ~/DATA/70deg_{$ver}_iter0.root -r --param ../params/invert_slope.ini -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_70deg_{$ver}.root -o ~/DATA/70deg_{$ver}_iter1.root -r --param ../params/invert_slope.ini -b

./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_60deg_{$ver}.root -o ~/DATA/60deg_c_{$ver}_iter0.root -r --param ../params/2by1.ini -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_60deg_{$ver}.root -o ~/DATA/60deg_c_{$ver}_iter1.root -r --param ../params/2by1.ini -b
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_60deg_{$ver}.root -o ~/DATA/60deg_d_{$ver}_iter0.root -r --param ../params/diag.ini -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_60deg_{$ver}.root- o ~/DATA/60deg_d_{$ver}_iter1.root -r --param ../params/diag.ini -b
~/DATA/60deg_c_{$ver}_iter1.root -r --param ../params/2by1.ini -b
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_60deg_{$ver}.root -o ~/DATA/60deg_{$ver}_iter0.root -r --param ../params/invert_slope.ini -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_60deg_{$ver}.root- o ~/DATA/60deg_{$ver}_iter1.root -r --param ../params/invert_slope.ini -b

./SpatialResol.exe -b -t 0 -i ~/DATA/tree_e-_4000MeV_50deg_{$ver}.root -o ~/DATA/50deg_d_{$ver}_iter0.root --param ../params/diag.ini -r; ./SpatialResol.exe -b -t 1 -i ~/DATA/tree_e-_4000MeV_50deg_{$ver}.root -o ~/DATA/50deg_d_{$ver}_iter1.root --param ../params/diag.ini -r
./SpatialResol.exe -b -t 0 -i ~/DATA/tree_e-_4000MeV_50deg_{$ver}.root -o ~/DATA/50deg_{$ver}_iter0.root -r; ./SpatialResol.exe -b -t 1 -i ~/DATA/tree_e-_4000MeV_50deg_{$ver}.root -o ~/DATA/50deg_{$ver}_iter1.root -r

./SpatialResol.exe -b -t 0 -i ~/DATA/tree_e-_4000MeV_45deg_{$ver}.root -o ~/DATA/45deg_d_{$ver}_iter0.root --param ../params/diag.ini -r; ./SpatialResol.exe -b -t 1 -i ~/DATA/tree_e-_4000MeV_45deg_{$ver}.root -o ~/DATA/45deg_d_{$ver}_iter1.root --param ../params/diag.ini -r
./SpatialResol.exe -b -t 0 -i ~/DATA/tree_e-_4000MeV_45deg_{$ver}.root -o ~/DATA/45deg_{$ver}_iter0.root -r; ./SpatialResol.exe -b -t 1 -i ~/DATA/tree_e-_4000MeV_45deg_{$ver}.root -o ~/DATA/45deg_{$ver}_iter1.root -r

./SpatialResol.exe -b -t 0 -i ~/DATA/tree_e-_4000MeV_40deg_{$ver}.root -o ~/DATA/40deg_d_{$ver}_iter0.root --param ../params/diag.ini -r; ./SpatialResol.exe -b -t 1 -i ~/DATA/40deg_{$ver}.root -o ~/DATA/40deg_d_{$ver}_iter1.root --param ../params/diag.ini -r
./SpatialResol.exe -b -t 0 -i ~/DATA/tree_e-_4000MeV_40deg_{$ver}.root -o ~/DATA/40deg_{$ver}_iter0.root -r; ./SpatialResol.exe -b -t 1 -i ~/DATA/tree_e-_4000MeV_40deg_{$ver}.root -o ~/DATA/40deg_{$ver}_iter1.root -r

./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_30deg_{$ver}.root -o ~/DATA/30deg_c_{$ver}_iter0.root -r --param ../params/2by1_invert.ini -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_30deg_{$ver}.root -o ~/DATA/30deg_c_{$ver}_iter1.root -r --param ../params/2by1_invert.ini -b
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_30deg_{$ver}.root -o ~/DATA/30deg_d_{$ver}_iter0.root -r --param ../params/diag.ini -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_30deg_{$ver}.root -o ~/DATA/30deg_d_{$ver}_iter1.root -r --param ../params/diag.ini -b
/SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_30deg_{$ver}.root -o ~/DATA/30deg_{$ver}_iter0.root -r -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_30deg_{$ver}.root -o ~/DATA/30deg_{$ver}_iter1.root -r -b

./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_20deg_{$ver}.root -o ~/DATA/20deg_c_{$ver}_iter0.root -r --param ../params/3by1_invert.ini -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_20deg_{$ver}.root -o ~/DATA/20deg_c_{$ver}_iter1.root -r --param ../params/3by1_invert.ini -b
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_20deg_{$ver}.root -o ~/DATA/20deg_d_{$ver}_iter0.root -r --param ../params/diag.ini -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_20deg_{$ver}.root -o ~/DATA/20deg_d_{$ver}_iter1.root -r --param ../params/diag.ini -b
./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_20deg_{$ver}.root -o ~/DATA/20deg_{$ver}_iter0.root -r -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_20deg_{$ver}.root -o ~/DATA/20deg_{$ver}_iter1.root -r -b

./SpatialResol.exe -t 0 -i ~/DATA/tree_e-_4000MeV_10deg_{$ver}.root -o ~/DATA/10deg_{$ver}_iter0.root -r -b; ./SpatialResol.exe -t 1 -i ~/DATA/tree_e-_4000MeV_10deg_{$ver}.root -o ~/DATA/10deg_{$ver}_iter1.root -r -b










