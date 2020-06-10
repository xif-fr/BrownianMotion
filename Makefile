CXX_FLAGS := -flto -O3 -std=c++11 -Wall
CXX_FLAGS_SFML := -DSFMLC01_WINDOW_UNIT=700 -DSFMLC01_WINDOW_HEIGHT=700 -DFONT_PATH=\"DejaVuSansMono.ttf\"
LD_FLAGS := -flto -lm -lfmt -fPIC
LD_FLAGS_SFML := -lsfml-graphics -lsfml-window -lsfml-system

brownian-motion-pysimul:
	echo "Compiling LJ-gas simulation core..."
	clang++ $(CXX_FLAGS) $(CXX_FLAGS_SFML) main.cpp pysimul-common.cpp pysimul-util.cpp vecs2.cpp $(LD_FLAGS) $(LD_FLAGS_SFML) -shared -o libpysimul.so
	echo "Building CFFI module..."
	python3 cffi-build.py
	rm *.o
	rm pysimul_ffi.c

brownian-motion-pysimul-headless: LD_FLAGS_SFML = 
brownian-motion-pysimul-headless: CXX_FLAGS_SFML = -DSIMUL_HEADLESS
brownian-motion-pysimul-headless: brownian-motion-pysimul

brownian-motion-standalone:
	clang++ $(CXX_FLAGS) -DMAIN_SOLO main.cpp pysimul-common.cpp pysimul-util.cpp vecs2.cpp $(LD_FLAGS) -o brownian-motion
	rm *.o

langevin-survival-pysimul:
	echo "Compiling Langevin survival simulation core..."
	clang++ $(CXX_FLAGS) -DSIMUL_HEADLESS langevin-survival.cpp pysimul-common.cpp pysimul-util.cpp vecs2.cpp $(LD_FLAGS) -shared -o libpysimul.so
	echo "Building CFFI module..."
	python3 cffi-build.py
	rm *.o
	rm pysimul_ffi.c
