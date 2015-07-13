#    This file is part of the WALT software
#
#    Copyright (C) 2015 University of Southern California
#                       Andrew D. Smith and Ting Chen
#
#    Authors: Haifeng Chen, Andrew D. Smith and Ting Chen
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(WALT)/src/smithlab_cpp
endif

app_subdirs=walt

all:
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(WALT) OPT=1; \
	done;

install:
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(WALT) OPT=1 install; \
	done;

test:
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(WALT) test; \
	done;

clean:
	@make -C smithlab_cpp OPT=1 clean
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(WALT) clean; \
	done;
