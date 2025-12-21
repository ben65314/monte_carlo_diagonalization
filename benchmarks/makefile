BINARY=unit_test

###DIRECTORIES###
LIBSDIR=../src/libs
INCDIRS=../src ../src/headers

###PARAMETERS###
CXX=g++ -std=c++20
OPT=-O3 -ftree-vectorize -lprofiler
# generate files that encode make rules for the .h dependencies
DEPFLAGS=-MP -MD
# automatically add the -I onto each include directory
WWD = -Wall -Wextra -Wpedantic -g
GCC = -llapack -lblas #-fopenmp #-lgfortran
CXXFLAGS=$(foreach D,$(INCDIRS),-I$(D)) $(OPT) $(DEPFLAGS) $(WWD)

# for-style iteration (foreach) and regular expression completions (wildcard)
CXXFILES=$(foreach D,$(LIBSDIR),$(wildcard $(D)/*.cpp))

# regular expression replacement
OBJECTS=$(patsubst %.cpp,%.o,$(CXXFILES))

DEPFILES=$(patsubst %.cpp,%.d,$(CXXFILES)) $(patsubst %.cpp,%.d,$(CXXFILESMAIN))

print-%  : ; @echo $* = $($*)
all: unit_test

unit_test: $(OBJECTS) ./unit_test.o
	$(CXX) -o $@ $^ $(GCC)

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(GCC)

OBJDIR:
	if [ ! -d $(OBJDIR) ]; then \
		mkdir $(OBJDIR); \
	fi

clean:
	rm -rf $(BINARY)
	rm -rf ./unit_test.o ./unit_test.d
retest:
	rm -rf unit_test
	rm -rf $(TESTDIR)/unit_test.o $(TESTDIR)/unit_test.d
	make

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
