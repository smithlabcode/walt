#
#    Copyright (C) 2015 University of Southern California
#                       Andrew D. Smith and Ting Chen
#
#    Authors: Haifeng Chen, Andrew D. Smith and Ting Chen
#
#    This file is part of the WALT.
#
#    WALT is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    WALT is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with WALT.  If not, see <http://www.gnu.org/licenses/>.
#

ifndef SRC_ROOT
SRC_ROOT=./../
endif

ifndef SMITHLAB_CPP
SMITHLAB_CPP=./../smithlab_cpp
endif

PROGS = makedb walt

CXX = g++
CFLAGS = -Wall -O3 -fPIC -fmessage-length=50 -fopenmp -D SEEDPATTERN3

TEST_DIR = $(SRC_ROOT)/test

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

SRCDB = reference.cpp
SRCMAP =  reference.cpp mapping.cpp paired.cpp
OBJDB = $(patsubst %.cpp,%.o,$(SRCDB))
OBJMAP = $(patsubst %.cpp,%.o,$(SRCMAP))

all:	$(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $<  -I$(SMITHLAB_CPP)

makedb : $(OBJDB) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o OptionParser.o)

walt : $(OBJMAP) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o OptionParser.o)

%: %.cpp
		$(CXX) $(CFLAGS) -o $@ $^ -I$(SMITHLAB_CPP)
install: all
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

test_%:	%
	@$(TEST_DIR)/$@ $(TEST_DIR)

test:	$(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
