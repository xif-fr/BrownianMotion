CXX_FLAGS := -flto -O3 -std=c++11 -Wall -Wno-switch
CXX_FLAGS_SFML := -DSFMLC01_WINDOW_UNIT=700 -DSFMLC01_WINDOW_HEIGHT=700 -DFONT_PATH=\"DejaVuSansMono.ttf\"
LD_FLAGS := -flto -lm -lfmt -fPIC
LD_FLAGS_SFML := -lsfml-graphics -lsfml-window -lsfml-system

build-cffi:
	echo "Building CFFI module..."
	python3 cffi-build.py
	rm *.o
	rm pysimul_ffi.c

brownian-motion-dynmol-pysimul: main.cpp
	echo "Compiling LJ-gas simulation core..."
	clang++ $(CXX_FLAGS) $(CXX_FLAGS_SFML) main.cpp pysimul-common.cpp pysimul-util.cpp vecs2.cpp $(LD_FLAGS) $(LD_FLAGS_SFML) -shared -o libpysimul.so
	make build-cffi

brownian-motion-dynmol-pysimul-headless: LD_FLAGS_SFML = 
brownian-motion-dynmol-pysimul-headless: CXX_FLAGS_SFML = -DSIMUL_HEADLESS
brownian-motion-dynmol-pysimul-headless: brownian-motion-dynmol-pysimul

brownian-motion-dynmol-standalone: main.cpp
	clang++ $(CXX_FLAGS) -DMAIN_SOLO main.cpp pysimul-common.cpp pysimul-util.cpp vecs2.cpp $(LD_FLAGS) -o brownian-motion
	rm *.o

brownian-motion-langevin-pysimul: langevin.cpp
	echo "Compiling Langevin brownian motion simulation core..."
	clang++ $(CXX_FLAGS) $(CXX_FLAGS_SFML) langevin.cpp pysimul-common.cpp pysimul-util.cpp vecs2.cpp $(LD_FLAGS) $(LD_FLAGS_SFML) -shared -o libpysimul.so
	make build-cffi

brownian-motion-langevin-pysimul-headless: LD_FLAGS_SFML = 
brownian-motion-langevin-pysimul-headless: CXX_FLAGS_SFML = -DSIMUL_HEADLESS
brownian-motion-langevin-pysimul-headless: brownian-motion-langevin-pysimul

langevin-survival-pysimul: langevin-survival.cpp
	echo "Compiling Langevin survival simulation core..."
	clang++ $(CXX_FLAGS) -DSIMUL_HEADLESS langevin-survival.cpp pysimul-common.cpp pysimul-util.cpp vecs2.cpp $(LD_FLAGS) -shared -o libpysimul.so
	make build-cffi

langevin-survival-simplified2d-pysimul: langevin-survival-simplified-2d.cpp
	echo "Compiling Langevin simplified 2D MFPT simulation core..."
	clang++ $(CXX_FLAGS) -DSIMUL_HEADLESS langevin-survival-simplified-2d.cpp pysimul-common.cpp pysimul-util.cpp vecs2.cpp $(LD_FLAGS) -shared -o libpysimul.so
	make build-cffi

langevin-generate: langevin-generate.cpp
	echo "Compiling standalone Langevin trajectory generator..."
	clang++ $(CXX_FLAGS) langevin-generate.cpp vecs2.cpp $(LD_FLAGS) -o langevin-generate

