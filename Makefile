# Source (C) Coire Cadeau 2007, all rights reserved.
# Modified and expanded by John Ming Ngo.

CC=g++
CCFLAGS=-std=c++17 -Wall -pedantic -O3 -g -fPIC
LDFLAGS=-lm

INCLUDE=include
SRC=src
COMMON=$(SRC)/common
SURFACE=$(SRC)/surface
SPECTRUM=$(SRC)/spectrum
SHAPE_SRC=shape_functions/src
SHAPE_LIB=shape_functions/lib
OBJDIR=build_objs

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    SHARED_EXT=.so
endif
ifeq ($(UNAME_S),Darwin)
    SHARED_EXT=.dylib
endif
ifeq ($(OS),Windows_NT)
    SHARED_EXT=.dll
    EXE_EXT=.exe
endif

PROGRAM=NS-SWORD$(EXE_EXT)

OBJ_NAMES=Main.o MainPrints.o OblModelBase.o PolyOblModelBase.o \
	OblDeflectionTOA.o MathFunctions.o Chi.o Hydrogen.o BlackBody.o \
	SpectrumManipulation.o SpectrumMapManipulation.o interp.o nrutil.o \
	Units.o matpack.o interp_functions.o utils.o
OBJ=$(addprefix $(OBJDIR)/,$(OBJ_NAMES))

SHAPE_FILES=$(wildcard $(SHAPE_SRC)/*.cpp)
SHARED_LIBS=$(patsubst $(SHAPE_SRC)/%.cpp,$(SHAPE_LIB)/%$(SHARED_EXT),$(SHAPE_FILES))
SHAPE_DEPS=$(OBJDIR)/OblModelBase.o $(OBJDIR)/MathFunctions.o $(OBJDIR)/Units.o \
	$(OBJDIR)/interp.o $(OBJDIR)/interp_functions.o $(OBJDIR)/nrutil.o

.PHONY: all main clean veryclean

all: $(PROGRAM)

main: $(PROGRAM)

$(PROGRAM): $(OBJ) $(SHARED_LIBS)
	$(CC) $(CCFLAGS) $(OBJ) $(LDFLAGS) -o $@

$(SHAPE_LIB)/%$(SHARED_EXT): $(SHAPE_SRC)/%.cpp $(SHAPE_DEPS) | $(SHAPE_LIB)
	$(CC) $(CCFLAGS) -shared $< $(OBJDIR)/OblModelBase.o $(OBJDIR)/MathFunctions.o \
		$(OBJDIR)/Units.o $(OBJDIR)/interp.o $(OBJDIR)/interp_functions.o \
		$(OBJDIR)/nrutil.o -o $@

$(OBJDIR)/Main.o: $(SRC)/main.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/MainPrints.o: $(SRC)/MainPrints.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/OblModelBase.o: $(SURFACE)/OblModelBase.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/PolyOblModelBase.o: $(SURFACE)/PolyOblModelBase.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/OblDeflectionTOA.o: $(SURFACE)/OblDeflectionTOA.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/MathFunctions.o: $(SURFACE)/MathFunctions.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/Chi.o: $(SURFACE)/Chi.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/Hydrogen.o: $(SPECTRUM)/Hydrogen.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/BlackBody.o: $(SPECTRUM)/BlackBody.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/SpectrumManipulation.o: $(SPECTRUM)/SpectrumManipulation.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/SpectrumMapManipulation.o: $(SPECTRUM)/SpectrumMapManipulation.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/interp.o: $(COMMON)/interp.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/nrutil.o: $(COMMON)/nrutil.c | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/Units.o: $(COMMON)/Units.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/matpack.o: $(COMMON)/matpack.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/interp_functions.o: $(COMMON)/interp_functions.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/utils.o: $(COMMON)/utils.cpp | $(OBJDIR)
	$(CC) $(CCFLAGS) -c $< -o $@

$(OBJDIR) $(SHAPE_LIB):
	mkdir -p $@

clean:
	rm -f $(OBJ)

veryclean: clean
	rm -f $(PROGRAM) $(SHARED_LIBS)
	rm -rf $(OBJDIR)/cache
