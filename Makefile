CXX ?= g++
CXXFLAGS ?= -std=c++20 -O3 -Wall -Wextra -Wpedantic -Wconversion -Iinclude

SRCDIR = src
INCDIR = include
BINDIR = bin

all: $(BINDIR)/cdbg

$(BINDIR):
	mkdir -p $(BINDIR)

$(BINDIR)/cdbg: $(SRCDIR)/main.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BINDIR)/cdbg_tests: tests/test_cdbg.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

test: $(BINDIR)/cdbg_tests
	./$(BINDIR)/cdbg_tests

clean:
	rm -rf $(BINDIR) build *.cdbg *.gfa

.PHONY: all test clean
