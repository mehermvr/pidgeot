export CC := /usr/local/bin/gcc
export CXX := /usr/local/bin/g++

.PHONY: build
build:
	cmake -G Ninja -S . -B build -DCMAKE_BUILD_TYPE=Debug -DCMAKE_VERBOSE_MAKEFILE=On
	cmake --build build

.PHONY: release_build
release_build:
	cmake -G Ninja -S . -B build -DCMAKE_BUILD_TYPE=Release
	cmake --build build

.PHONY: clean
clean:
	rm -rf build
