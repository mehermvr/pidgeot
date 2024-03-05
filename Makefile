export CC := /usr/bin/clang-18
export CXX := /usr/bin/clang++-18

.PHONY: build
build:
	cmake -G Ninja -S . -B build -DCMAKE_BUILD_TYPE=Debug
	cmake --build build -j $(shell nproc)

.PHONY: release_build
release_build:
	cmake -G Ninja -S . -B build -DCMAKE_BUILD_TYPE=Release
	cmake --build build -j $(shell nproc)

.PHONY: clean
clean:
	rm -rf build
