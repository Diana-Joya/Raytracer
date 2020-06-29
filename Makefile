
CLANG_FLAGS = -std=c++17 -Wall -O -g
PNG_FLAGS = `libpng-config --cflags --ldflags`
GTEST_FLAGS = -lpthread -lgtest_main -lgtest  -lpthread

default: fast

all: fast slow

headers: gfxalgebra.hpp gfximage.hpp gfxnumeric.hpp gfxpng.hpp gfxraytrace.hpp rayson.hpp

libraries: /usr/lib/libgtest.a /usr/include/png++/png.hpp

mrraytracer: headers libraries mrraytracer.cpp
	clang++ ${CLANG_FLAGS} ${PNG_FLAGS} mrraytracer.cpp -o mrraytracer

%.png : %.json mrraytracer
	./mrraytracer $< $@

fast: scene_gtri_ortho_flat.png \
	    scene_gtri_persp_flat.png \
			scene_gtri_ortho_phong.png \
			scene_gtri_persp_phong.png \
			scene_2spheres_ortho_flat.png \
			scene_2spheres_persp_flat.png \
			scene_2spheres_ortho_phong.png \
			scene_2spheres_persp_phong.png

slow: teatime.png

/usr/lib/libgtest.a:
	@echo -e "google test library not installed\n"
	@echo -e "Installing libgtest-dev. Please provide the password when asked\n"
	@sudo apt-get -y install libgtest-dev cmake
	@sudo apt-get install cmake # install cmake
	@echo -e "\nConfiguring libgtest-dev\n"
	@cd /usr/src/gtest; sudo cmake CMakeLists.txt; sudo make; sudo cp *.a /usr/lib
	@echo -e "Finished installing google test library\n"

/usr/include/png++/png.hpp:
	@echo -e "libpng++ library not installed\n"
	@echo -e "Installing libpng++-dev. Please provide the password when asked\n"
	@sudo apt-get -y install libpng++-dev

clean:
		rm -f mrraytracer *png
