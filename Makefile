# programs
TARGET := lcpan
SRCS := $(wildcard *.c)

OBJS := $(SRCS:.c=.o)
DOBJS := $(SRCS:.c=.dbg.o)

# directories
CURRENT_DIR := $(shell pwd)

GXX ?= gcc
OPT_FLAGS=-O3 -march=native -ftree-vectorize
DBG_OPT_FLAGS=-O0 -g   -fsanitize=address -fsanitize=undefined 
CXXFLAGS =  -Wall -Wextra -Wpedantic  -lz -Wno-interference-size  
ifeq (${GXX},g++)
CXXFLAGS := ${CXXFLAGS} -std=c++23
endif

# object files that need lcptools
#LCPTOOLS_CXXFLAGS := -I$(CURRENT_DIR)/lcptools/include
LCPTOOLS_LDFLAGS :=  -Wl,-rpath,$(CURRENT_DIR)/lcptools/lib -lz

all: $(TARGET) $(TARGET)_dbg


lcptools_ho.h: lcptools/lcptools_ho.h
	cp lcptools/lcptools_ho.h .

$(TARGET)_dbg: $(DOBJS)
	$(GXX) $(DBG_OPT_FLAGS) $(CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -o $@ $^ $(LCPTOOLS_LDFLAGS) -lm -lpthread
	rm $(DOBJS)
	@mkdir -p bin
	mv $(TARGET)_dbg bin

$(TARGET): $(OBJS)
	$(GXX) $(OPT_FLAGS) $(CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -o $@ $^ $(LCPTOOLS_LDFLAGS) -lm -lpthread
	rm $(OBJS)
	@mkdir -p bin
	mv $(TARGET) bin

%.o: %.c lcptools_ho.h 
	$(GXX) $(OPT_FLAGS) $(CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -c $< -o $@

%.dbg.o: %.c
	$(GXX) $(DBG_OPT_FLAGS) $(CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -c $< -o $@
install:
	@echo "Installing lcptools"
	cd lcptools && make PREFIX=. install

clean:
	rm -f $(TARGET) $(OBJS) $(TARGET)_dbg $(DOBJS)
