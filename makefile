BINARY=mcd_solver subspace_finder qMatrixGreenCompute graphGreen av_combine

###DIRECTORIES###
LIBSDIR=src/libs
MAINDIR = src/mainCode
INCDIRS=src src/headers

###PARAMETERS###
CXX=g++ -std=c++20
OPT=-O3 -ftree-vectorize -lprofiler
# generate files that encode make rules for the .h dependencies
DEPFLAGS=-MP -MD
# automatically add the -I onto each include directory
WWD = -Wall -Wextra -Wpedantic -g -Wno-unused-parameter 
GCC = -llapack -llapacke -lblas -fopenmp #-lgfortran 
CXXFLAGS=-I$(foreach D,$(INCDIRS),-I$(D)) $(OPT) $(DEPFLAGS) $(WWD)  

# for-style iteration (foreach) and regular expression completions (wildcard)
CXXFILES=$(foreach D,$(LIBSDIR),$(wildcard $(D)/*.cpp))
CXXFILESMAIN=$(foreach D,$(MAINDIR),$(wildcard $(D)/*.cpp))

# regular expression replacement
OBJECTS=$(patsubst %.cpp,%.o,$(CXXFILES))

MAIN_OBJECTS = $(patsubst %.cpp,%.o,$(CXXFILESMAIN))

DEPFILES=$(patsubst %.cpp,%.d,$(CXXFILES)) $(patsubst %.cpp,%.d,$(CXXFILESMAIN))

GREENGRAPHFILE=$(MAINDIR)

print-%  : ; @echo $* = $($*)

all: build mcd_solver subspace_finder qMatrixGreenCompute graphGreen av_combine
	mv mcd_solver build/
	mv subspace_finder build/

mcd_solver: $(OBJECTS) $(MAINDIR)/mcd_solver.o
	$(CXX) -o $@ $^ $(GCC)
subspace_finder: $(OBJECTS) $(MAINDIR)/subspace_finder.o
	$(CXX) -o $@ $^ $(GCC)
	
qMatrixGreenCompute:
	cp -p $(GREENGRAPHFILE)/qMatrixGreenCompute.py ./build/qMatrixGreenCompute.py
graphGreen:
	cp -p $(GREENGRAPHFILE)/graphGreen.py ./build/graphGreen.py
av_combine:
	cp -p $(GREENGRAPHFILE)/av_combine.py ./build/av_combine.py

# only want the .c file dependency here, thus $< instead of $^.
#
%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(GCC)
	
OBJDIR:
	if [ ! -d $(OBJDIR) ]; then \
		mkdir $(OBJDIR); \
	fi

build:
	rm -rf build
	mkdir -p build

clean:
	rm -rf $(BINARY) $(OBJECTS) $(DEPFILES) $(MAIN_OBJECTS) $(OBJDIR)
	rm -rf build



# shell commands are a set of keystrokes away
distribute: clean
	tar zcvf dist.tgz *

# @ silences the printing of the command
# $(info ...) prints output
diff:
	$(info The status of the repository, and the volume of per-file changes:)
	@git status
	@git diff --stat

# include the dependencies
-include $(DEPFILES)

# add .PHONY so that the non-targetfile - rules work even if a file with the same name exists.
.PHONY: all clean distribute diff
