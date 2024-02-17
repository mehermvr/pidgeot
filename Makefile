.PHONY: build
build:
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
	cmake --build build -j $(shell nproc)

.PHONY: release_build
release_build:
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
	cmake --build build -j $(shell nproc)

.PHONY: clean
clean:
	rm -rf build
