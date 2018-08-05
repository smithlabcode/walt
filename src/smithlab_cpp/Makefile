# Makefile from smithlab_cpp C++ code library
#
# Copyright (C) 2010-2018 University of Southern California and
#                         Andrew D. Smith
#
# Authors: Andrew D. Smith
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))

CXX = g++
CXXFLAGS = -Wall -std=c++11
OPTFLAGS = -O2
DEBUGFLAGS = -g

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(OBJECTS)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

clean:
	@-rm -f *.o *~
.PHONY: clean
