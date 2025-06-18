# programs
TARGET := lcpan
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)

# directories
CURRENT_DIR := $(shell pwd)

GXX ?= gcc
CXXFLAGS = -O3 -g -Wall -Wextra -Wpedantic -march=native -ftree-vectorize -lz 
ifeq (${GXX},g++)
CXXFLAGS := ${CXXFLAGS} -std=c++23
endif

# object files that need lcptools
LCPTOOLS_CXXFLAGS := -I$(CURRENT_DIR)/lcptools/include
LCPTOOLS_LDFLAGS := -L$(CURRENT_DIR)/lcptools/lib -llcptools -Wl,-rpath,$(CURRENT_DIR)/lcptools/lib -lz

$(TARGET): $(OBJS)
	$(GXX) $(CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -o $@ $^ $(LCPTOOLS_LDFLAGS) -lm -lpthread
	rm $(OBJS)
	@mkdir -p bin
	mv $(TARGET) bin

%.o: %.c
	$(GXX) $(CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -c $< -o $@

install:
	@echo "Installing lcptools"
	cd lcptools && make PREFIX=. install

clean:
	rm -f $(TARGET) $(OBJS)
