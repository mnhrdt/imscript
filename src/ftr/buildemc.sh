BIN=miniev
BIN=jmgs
emcc -c $BIN.c
emcc -c ftr_emscripten.c
emcc $BIN.o ftr_emscripten.o -o $BIN.html -s SINGLE_FILE=1
