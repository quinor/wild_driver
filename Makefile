all: wild_driver

.PHONY: all clean

wild_driver: src/main.cc src/external/CLI11.hpp src/external/linalg.h src/external/nanosvg.h Makefile
	g++ -Wall -Wextra -std=c++11 src/main.cc -o wild_driver

clean:
	rm wild_driver
