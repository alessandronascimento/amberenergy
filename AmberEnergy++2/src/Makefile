#!/bin/bash

CC=g++
CFLAGS=-c -O3 -ffast-math -w -g3 -I../lib/zlib-1.2.11 -I../lib/gzstream-master -I../lib/netcdf/include
LDFLAGS= -lgzstream -lz -lnetcdf -lnetcdf_c++ -L../lib/netcdf/lib -L../lib/zlib-1.2.11 -L../lib/gzstream-master 
SOURCES=COORD.cpp  ENERGY.cpp  main.cpp  PRMTOP.cpp  RunEngine.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=AmberEnergy++

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(OBJECTS)
