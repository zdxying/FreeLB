# ============================================================================
# Makefile for xcore MemoryPool library
#
# Builds a static library (libxcore_mem.a) and/or a shared library
# (libxcore_mem.so) that other programs can link against.
#
# Usage:
#   make            build both static and shared libraries (release)
#   make static     build static library only
#   make shared     build shared library only
#   make DEBUG=1    debug build (-O0 -g, defines DEBUG and MEMPOOL_TRACK)
#   make clean      remove build artifacts
#   make install    install headers + libraries into $(PREFIX) (default /usr/local)
#   make uninstall  remove installed files
#
# Release and debug builds live in separate directories (build/release and
# build/debug), so switching `make` <-> `make DEBUG=1` never reuses a stale
# object file compiled with different flags.
# ============================================================================

CXX      ?= g++
AR       ?= ar

LIB_NAME := xcore
PREFIX   ?= /usr/local

# build configuration ---------------------------------------------------------

ifeq ($(DEBUG),1)
CONFIG := debug
CXXFLAGS = -O0 -g -std=c++17 -Wall -Wextra -fPIC
CPPFLAGS = -Isrc -DDEBUG -DMEMPOOL_TRACK
else
CONFIG := install
CXXFLAGS = -O2 -std=c++17 -march=native -Wall -Wextra -fPIC
CPPFLAGS = -Isrc
endif

BUILD_DIR := build/$(CONFIG)
OBJ_DIR := $(BUILD_DIR)/obj
LIB_DIR := $(BUILD_DIR)/lib

SRC_DIR := src
SRCS := $(shell find $(SRC_DIR) -name '*.cpp')
OBJS := $(patsubst src/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

STATIC_LIB := $(LIB_DIR)/lib$(LIB_NAME).a
SHARED_LIB := $(LIB_DIR)/lib$(LIB_NAME).so

.PHONY: all static shared clean install uninstall

all: static shared

static: $(STATIC_LIB)

shared: $(SHARED_LIB)

$(STATIC_LIB): $(OBJS)
	@mkdir -p $(LIB_DIR)
	$(AR) rcs $@ $^

$(SHARED_LIB): $(OBJS)
	@mkdir -p $(LIB_DIR)
	$(CXX) -shared -o $@ $^

$(OBJ_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rf build

install: all
	install -d $(PREFIX)/include/xcore/memory $(PREFIX)/include/xcore/util $(PREFIX)/lib
	install -m 644 src/xcore/memory/memoryPool.h $(PREFIX)/include/xcore/memory/
	install -m 644 src/xcore/memory/poolAllocator.h $(PREFIX)/include/xcore/memory/
	install -m 644 src/xcore/util/assert.h $(PREFIX)/include/xcore/util/
	install -m 644 $(STATIC_LIB) $(PREFIX)/lib/
	install -m 755 $(SHARED_LIB) $(PREFIX)/lib/

uninstall:
	rm -f $(PREFIX)/include/xcore/memory/memoryPool.h
	rm -f $(PREFIX)/include/xcore/memory/poolAllocator.h
	rm -f $(PREFIX)/include/xcore/util/assert.h
	rm -f $(PREFIX)/lib/lib$(LIB_NAME).a
	rm -f $(PREFIX)/lib/lib$(LIB_NAME).so

-include $(DEPS)