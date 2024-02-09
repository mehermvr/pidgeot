.PHONY: build
build:
	cmake -S . -B build
	cmake --build build -j $(nproc)

.PHONY: clean
clean:
	rm -rf build
