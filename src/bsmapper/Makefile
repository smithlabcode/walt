ifndef SRC_ROOT
SRC_ROOT=./../
endif

ifndef SMITHLAB_CPP
SMITHLAB_CPP=./../smithlab_cpp
endif

PROGS = makedb bsmapper

CXX = g++
CFLAGS = -Wall -O3 -fPIC -fmessage-length=50 -pg

TEST_DIR = $(SRC_ROOT)/test

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

SRCDB = reference.cpp
SRCMAP =  reference.cpp mapping.cpp
OBJDB = $(patsubst %.cpp,%.o,$(SRCDB))
OBJMAP = $(patsubst %.cpp,%.o,$(SRCMAP))

all:	$(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $<  -I$(SMITHLAB_CPP)

makedb : $(OBJDB) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o OptionParser.o)

bsmapper : $(OBJMAP) \
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
