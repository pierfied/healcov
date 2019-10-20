all: healcov

.PHONY: healcov
healcov:
	mkdir -p cmake_build && \
	cmake -Bcmake_build -H. && \
	cd cmake_build && \
	make && \
	cp healcov ../healcov/healcov

.PHONY: clean
clean:
	rm -rf cmake_build
	rm -rf healcov/healcov