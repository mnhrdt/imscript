emcc -c miniev.c
emcc -c ftr_emscripten.c
emcc miniev.o ftr_emscripten.o -o miniev.html -s SINGLE_FILE=1
