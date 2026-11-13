# ChemLarp — make-based build for macOS, Linux and Windows (MSYS2/MinGW-w64).
# On Windows with Visual Studio, use CMake instead (scripts/build.ps1).
CXX      ?= clang++
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wshadow
CPPFLAGS += -Iinclude

ifeq ($(OS),Windows_NT)
  # MinGW-w64 (MSYS2): self-contained DLL and executables, no libstdc++ DLLs needed.
  CXX         := g++
  EXE         := .exe
  SHARED      := build/chemlarp.dll
  SHARED_FLAG := -shared -Wl,--out-implib,build/libchemlarp.dll.a
  LDFLAGS     += -static-libgcc -static-libstdc++ -static
else
  UNAME    := $(shell uname -s)
  CXXFLAGS += -fPIC
  EXE      :=
  ifeq ($(UNAME),Darwin)
    SHARED      := build/libchemlarp.dylib
    SHARED_FLAG := -dynamiclib -install_name @rpath/libchemlarp.dylib
  else
    SHARED      := build/libchemlarp.so
    SHARED_FLAG := -shared
  endif
endif

LIB_SRC  := $(wildcard src/*.cpp)
LIB_OBJ  := $(LIB_SRC:src/%.cpp=build/obj/%.o)
TEST_SRC := $(wildcard tests/*.cpp)
TEST_OBJ := $(TEST_SRC:tests/%.cpp=build/obj/tests/%.o)
CLI      := build/chemlarp$(EXE)
TESTS    := build/chemlarp_tests$(EXE)

.PHONY: all test clean

all: build/libchemlarp.a $(SHARED) $(CLI) $(TESTS)

build/obj/%.o: src/%.cpp $(wildcard include/chemlarp/*.hpp include/chemlarp/*.h src/*.hpp src/*.inc)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

build/obj/tests/%.o: tests/%.cpp $(wildcard include/chemlarp/*.hpp tests/*.hpp)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

build/obj/cli/main.o: src/cli/main.cpp $(wildcard include/chemlarp/*.hpp)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

build/libchemlarp.a: $(LIB_OBJ)
	ar rcs $@ $^

$(SHARED): $(LIB_OBJ)
	$(CXX) $(SHARED_FLAG) $(LDFLAGS) $^ -o $@

$(CLI): build/obj/cli/main.o build/libchemlarp.a
	$(CXX) $(LDFLAGS) $^ -o $@

$(TESTS): $(TEST_OBJ) build/libchemlarp.a
	$(CXX) $(LDFLAGS) $^ -o $@

test: $(TESTS)
	./$(TESTS)

clean:
	rm -rf build
