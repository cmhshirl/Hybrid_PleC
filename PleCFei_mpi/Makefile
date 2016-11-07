#the compiler
CC = mpic++ 

# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbounds-check
FCFLAGS = -O3 
# flags forall (e.g. look for system .mod files, required in gfortran)
# FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = Stochastic 
SOURCES = Stochastic.cpp PleC.cpp HillReaction.cpp Reaction.cpp Species.cpp \
          CatReaction.cpp CmbReaction.cpp DecmpstReaction.cpp \
          LftDiffusion.cpp RgtDiffusion.cpp

OBJECTS = $(SOURCES:.cpp=.o)

# "make" builds all
all: $(PROGRAMS)

Stochastic : Stochastic.o PleC.o HillReaction.o Reaction.o Species.o
Stochastic.o : PleC.h 
PleC.o : Reaction.h CatReaction.h CmbReaction.h DecmpstReaction.h \
         HillReaction.h LftDiffusion.h RgtDiffusion.h
Reaction.o : Species.h
CatReaction.o : Reaction.h Species.h
CmbReaction.o : Reaction.h Species.h
DecmpstReaction.o : Reaction.h Species.h
HillReaction.o : Reaction.h Species.h
LftDiffusion.o : Reaction.h Species.h
RgtDiffusion.o : Reaction.h Species.h

$(PROGRAMS): $(OBJECTS)
	$(CC) $(FCFLAGS) $(OBJECTS) -o $@

%.o: %.cpp
	$(CC) $(FCFLAGS) -c $<


.PHONY : clean
clean : 
	rm -f *.o *~ 

veryclean: clean
	rm -f $(PROGRAMS) 


