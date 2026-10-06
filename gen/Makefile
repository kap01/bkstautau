PREFIX   ?= $(HOME)/bkstautau/install
CXX      ?= g++
FJCONFIG := $(PREFIX)/bin/fastjet-config
CXXFLAGS := -O2 -std=c++17 -pthread -I$(PREFIX)/include $(shell $(FJCONFIG) --cxxflags)
LDFLAGS  := -L$(PREFIX)/lib -Wl,-rpath,$(PREFIX)/lib -lpythia8 -ldl $(shell $(FJCONFIG) --libs)

all: gen_bkstautau

gen_bkstautau: gen_bkstautau.cc
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

run: gen_bkstautau
	./gen_bkstautau bkstautau.cmnd bkstautau.csv

clean:
	rm -f gen_bkstautau

.PHONY: all run clean
