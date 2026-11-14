# ChemLarp — make-based build for macOS, Linux and Windows (MSYS2/MinGW-w64).
# On Windows with Visual Studio, use CMake instead (scripts/build.ps1).
# BUILD: output directory (a separate one for builds with other flags, e.g. the Python wheel).
BUILD    ?= build
CXX      ?= clang++
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wshadow
CPPFLAGS += -Iinclude

ifeq ($(OS),Windows_NT)
  # MinGW-w64 (MSYS2): self-contained DLL and executables, no libstdc++ DLLs needed.
  CXX         := g++
  EXE         := .exe
  SHARED      := $(BUILD)/chemlarp.dll
  SHARED_FLAG := -shared -Wl,--out-implib,$(BUILD)/libchemlarp.dll.a
  LDFLAGS     += -static-libgcc -static-libstdc++ -static
else
  UNAME    := $(shell uname -s)
  CXXFLAGS += -fPIC
  EXE      :=
  ifeq ($(UNAME),Darwin)
    SHARED      := $(BUILD)/libchemlarp.dylib
    SHARED_FLAG := -dynamiclib -install_name @rpath/libchemlarp.dylib
  else
    SHARED      := $(BUILD)/libchemlarp.so
    SHARED_FLAG := -shared
    CXXFLAGS    += -pthread
    LDFLAGS     += -pthread
  endif
endif

# ARCH_FLAGS (e.g. -arch x86_64 when cross-compiling on macOS) and EXTRA_LDFLAGS (e.g. -static-libstdc++) are
# added to the platform's flags; the Python wheel build uses them.
CXXFLAGS += $(ARCH_FLAGS)
LDFLAGS  += $(ARCH_FLAGS) $(EXTRA_LDFLAGS)

LIB_SRC  := $(wildcard src/*.cpp)
LIB_OBJ  := $(LIB_SRC:src/%.cpp=$(BUILD)/obj/%.o)
TEST_SRC := $(wildcard tests/*.cpp)
TEST_OBJ := $(TEST_SRC:tests/%.cpp=$(BUILD)/obj/tests/%.o)
CLI      := $(BUILD)/chemlarp$(EXE)
TESTS    := $(BUILD)/chemlarp_tests$(EXE)

.PHONY: all test validate install clean

all: $(BUILD)/libchemlarp.a $(SHARED) $(CLI) $(TESTS)

$(BUILD)/obj/%.o: src/%.cpp $(wildcard include/chemlarp/*.hpp include/chemlarp/*.h src/*.hpp src/*.inc)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# The tests link the static library, so the C API is declared without dllimport (as in CMake).
$(BUILD)/obj/tests/%.o: tests/%.cpp $(wildcard include/chemlarp/*.hpp include/chemlarp/*.h tests/*.hpp)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) -DCHEMLARP_STATIC $(CXXFLAGS) -c $< -o $@

$(BUILD)/obj/cli/main.o: src/cli/main.cpp $(wildcard include/chemlarp/*.hpp)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD)/libchemlarp.a: $(LIB_OBJ)
	ar rcs $@ $^

$(SHARED): $(LIB_OBJ)
	$(CXX) $(SHARED_FLAG) $(LDFLAGS) $^ -o $@

$(CLI): $(BUILD)/obj/cli/main.o $(BUILD)/libchemlarp.a
	$(CXX) $(LDFLAGS) $^ -o $@

$(TESTS): $(TEST_OBJ) $(BUILD)/libchemlarp.a
	$(CXX) $(LDFLAGS) $^ -o $@

test: $(TESTS)
	./$(TESTS)

# Every validation suite that can run here, and a report (see scripts/validate.py for the data layout):
#   make validate DATA=path/to/reference-data PY=path/to/venv/bin/python
validate: all
	python3 scripts/validate.py $(if $(DATA),--data $(DATA),) $(if $(PY),--python $(PY),) \
	  --report $(BUILD)/validation-report.md --json $(BUILD)/validation-report.json

# make install [PREFIX=/usr/local] [DESTDIR=...]: the CLI, the shared library, the C/C++ headers and the
# data files (example mechanisms, reference structures) under PREFIX.
PREFIX ?= /usr/local
install: $(SHARED) $(CLI)
	install -d $(DESTDIR)$(PREFIX)/bin $(DESTDIR)$(PREFIX)/lib $(DESTDIR)$(PREFIX)/include/chemlarp $(DESTDIR)$(PREFIX)/share/chemlarp
	install -m 755 $(CLI) $(DESTDIR)$(PREFIX)/bin/
	install -m 644 $(SHARED) $(DESTDIR)$(PREFIX)/lib/
	install -m 644 include/chemlarp/*.h include/chemlarp/*.hpp $(DESTDIR)$(PREFIX)/include/chemlarp/
	cp -R data/. $(DESTDIR)$(PREFIX)/share/chemlarp/

clean:
	rm -rf $(BUILD)
