all: wild_driver

.PHONY: all clean

EXT_HEADERS := src/external/CLI11.hpp src/external/linalg.h src/external/nanosvg.h src/external/json.hpp

HEADERS := src/svg.hh src/math.hh src/json_output.hh
SOURCES := src/main.cc src/svg.cc src/nanosvg_impl.cc src/json_output.cc

wild_driver: $(SOURCES) $(HEADERS) $(EXT_HEADERS) Makefile
	g++ -Wall -Wextra -std=c++14 $(SOURCES) -o wild_driver

clean:
	rm wild_driver
