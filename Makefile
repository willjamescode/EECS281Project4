# A Makefile that specifies a couple different ways to compile

EXECUTABLE = tsp
MAINFILE = tspSolver.cpp

# designate which compiler to use
COMPILER = g++

#Compiler flags used by default
COMPFLAGS = -std=c++11 

# list of cpp source files
SOURCES = $(wildcard *.cpp)

# list of object files
OBJ = $(SOURCES:%.cpp=%.o)

# make release - compiles "all" with $(COMPFLAGS) and -O3 flag
# also specifies NDEBUG, so asserts will not check
release: COMPFLAGS += -O3
release: all

# make profile - compiles "all" with $(COMPFLAGS) and -pg flag
profile: COMPFLAGS += -pg
profile: clean all

# make debug - compiles "all" with $(COMPFLAGS) and -g flag
# also specifies DEBUG so "#ifdef DEBUG /*...*/ #endif" works
debug: COMPFLAGS += -g3 -DDEBUG
debug: clean all

# strings all the files together into an executable
all: $(OBJ)
ifeq ($(EXECUTABLE), executable)
	$(COMPILER) $(COMPFLAGS) $(OBJ)
else
	$(COMPILER) $(COMPFLAGS) $(OBJ) -o $(EXECUTABLE)
endif

# specifies how to create objects 
%.o:
	$(COMPILER) $(COMPFLAGS) -c $*.cpp

# make clean - removes all objects and executable
clean:
	rm -f $(OBJ) $(EXECUTABLE)

# disables built in rules
.SUFFIXES:
# target that doesn't create any files
.PHONY: all release debug profile clean
