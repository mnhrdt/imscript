export CFLAGS=-march=native -Os

make -C src -j 3
make -C src_ftr -j 3
echo cp src/+([a-z]) src_ftr/+([a-z]) bin
