# Makefile

BINDIR := bin
LIBDIR := lib

CCFLAGS :=  -Wextra -Wfloat-equal -Wundef -Werror -fverbose-asm  -Wshadow -Wpointer-arith -Wcast-align -Wconversion
DEBUG := -g -pedantic -Wall
CC := mpic++ -std=c++14 -O2 -pedantic


VPATH = main:src

SRC := $(wildcard src/*.cc)
OBJ := $(patsubst %.cc, $(BINDIR)/%.o, $(notdir $(SRC)))
INC := $(wildcard src/*.hh)


mpiAd: MPIAdvance.exe
	@cd bin; echo 'running program...\n \nOutput Results:'; mpirun ./MPIAdvance.exe;

Advance: Advancetest.exe
	@cd bin; echo 'running program...\n \nOutput Results:'; ./Advancetest.exe;

test: test.exe
	@cd bin; echo "running program... \n Output Results:" ; mpirun ./test.exe;

Ctest: Ctest.exe
	@cd bin; echo 'running program...\nOutput Results:'; mpirun -np 4 ./Ctest.exe;

lib: $(LIBDIR)/libPIC.a

$(LIBDIR)/libPIC.a: $(OBJ) 
	@echo make lib...
	ar rv $@ $^
	ranlib $@

%.exe: $(BINDIR)/%.o $(LIBDIR)/libPIC.a 
	@echo compilink and linking... 
	$(CC) -I src $< -o $(BINDIR)/$@ -L lib -l PIC

$(BINDIR)/%.o: %.cc | $(INC)
	@echo compiling... $<
	$(CC) -I src -c $< -o $@

######### clean

tilde := $(wildcard */*~) $(wildcard *~)
exe := $(wildcard */*.exe) $(wildcard *.exe)
obj := $(wildcard */*.o) $(wildcard *.o)  $(wildcard */*.pcm) $(wildcard */*.d)
mylibs := $(wildcard */*.so) $(wildcard */*.a)

txtOut := $(wildcard results/*.txt)

clean:
	@echo cleaning dir...
	rm -f $(exe) $(obj) $(tilde) $(mylibs)

OutClean:
	@echo cleaning results...
	@rm -f $(txtOut)$