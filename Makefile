all: healcov

.PHONY: healcov
healcov:
	mkdir -p cmake_build && \
	cmake -Bcmake_build -H. && \
	cd cmake_build && \
	make && \
	cp libhealcov.so ../healcov/libhealcov.so

.PHONY: clean
clean:
	rm -rf cmake_build
	rm -rf healcov/libhealcov.so