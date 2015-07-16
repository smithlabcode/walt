#
#    Copyright (C) 2015 University of Southern California
#                       Andrew D. Smith and Ting Chen
#
#    Authors: Haifeng Chen, Andrew D. Smith and Ting Chen
#
#    This file is part of the WALT
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

WALT = $(shell pwd)

all:
	@make -C src WALT=$(WALT) OPT=1

install:
	@make -C src WALT=$(WALT) OPT=1 install

test:
	@make -C src WALT=$(WALT) test
.PHONY: test

clean:
	@make -C src WALT=$(WALT) clean
.PHONY: clean

distclean: clean
	@rm -rf $(WALT)/bin
	@rm -rf $(WALT)/include
.PHONY: distclean
