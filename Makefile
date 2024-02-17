.PHONY: build
build:
	cmake -S . -B build
	cmake --build build -j $(shell nproc)

.PHONY: clean
clean:
	rm -rf build
