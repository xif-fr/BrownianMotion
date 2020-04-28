CXX_FLAGS := -flto -O3 -std=c++11 -Wall -DSFMLC01_WINDOW_UNIT=700 -DSFMLC01_WINDOW_HEIGHT=700 -DFONT_PATH=\"DejaVuSansMono.ttf\"
LD_FLAGS := -flto -lsfml-graphics -lsfml-window -lsfml-system -lm -lfmt -fPIC

brownian-motion-pysimul:
	echo "Compiling simulation core..."
	clang++ $(CXX_FLAGS) main.cpp pysimul-common.cpp vecs2.cpp $(LD_FLAGS) -shared -o libpysimul.so
	echo "Building CFFI module..."
	python3 cffi-build.py
	rm *.o
	rm pysimul_ffi.c

brownian-motion-standalone:
	clang++ $(CXX_FLAGS) -DMAIN_SOLO main.cpp pysimul-common.cpp vecs2.cpp $(LD_FLAGS) -o brownian-motion
	rm *.o
