# OmniChem Core — make-based build for macOS, Linux and Windows (MSYS2/MinGW-w64).
# On Windows with Visual Studio, use CMake instead (scripts/build.ps1).
CXX      ?= clang++
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wshadow
CPPFLAGS += -Iinclude

ifeq ($(OS),Windows_NT)
  # MinGW-w64 (MSYS2): self-contained DLL and executables, no libstdc++ DLLs needed.
  CXX         := g++
  EXE         := .exe
  SHARED      := build/omnichem.dll
  SHARED_FLAG := -shared -Wl,--out-implib,build/libomnichem.dll.a
  LDFLAGS     += -static-libgcc -static-libstdc++ -static
else
  UNAME    := $(shell uname -s)
  CXXFLAGS += -fPIC
  EXE      :=
  ifeq ($(UNAME),Darwin)
    SHARED      := build/libomnichem.dylib
    SHARED_FLAG := -dynamiclib -install_name @rpath/libomnichem.dylib
  else
    SHARED      := build/libomnichem.so
    SHARED_FLAG := -shared
  endif
endif

LIB_SRC  := $(wildcard src/*.cpp)
LIB_OBJ  := $(LIB_SRC:src/%.cpp=build/obj/%.o)
TEST_SRC := $(wildcard tests/*.cpp)
TEST_OBJ := $(TEST_SRC:tests/%.cpp=build/obj/tests/%.o)
CLI      := build/omnichem$(EXE)
TESTS    := build/omnichem_tests$(EXE)

.PHONY: all test clean

all: build/libomnichem.a $(SHARED) $(CLI) $(TESTS)

build/obj/%.o: src/%.cpp $(wildcard include/omnichem/*.hpp include/omnichem/*.h)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

build/obj/tests/%.o: tests/%.cpp $(wildcard include/omnichem/*.hpp tests/*.hpp)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

build/obj/cli/main.o: src/cli/main.cpp $(wildcard include/omnichem/*.hpp)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

build/libomnichem.a: $(LIB_OBJ)
	ar rcs $@ $^

$(SHARED): $(LIB_OBJ)
	$(CXX) $(SHARED_FLAG) $(LDFLAGS) $^ -o $@

$(CLI): build/obj/cli/main.o build/libomnichem.a
	$(CXX) $(LDFLAGS) $^ -o $@

$(TESTS): $(TEST_OBJ) build/libomnichem.a
	$(CXX) $(LDFLAGS) $^ -o $@

test: $(TESTS)
	./$(TESTS)

clean:
	rm -rf build
