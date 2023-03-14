# Makefile

BINDIR := bin
LIBDIR := lib

CCFLAGS :=  -Wextra -Wfloat-equal -Wundef -Werror -fverbose-asm  -Wshadow -Wpointer-arith -Wcast-align -Wconversion
DEBUG := -g -pedantic -Wall
CC := mpic++ -std=c++14 -O2 -pedantic


VPATH = main:src

SRC := $(wildcard src/*.C)
OBJ := $(patsubst %.C, $(BINDIR)/%.o, $(notdir $(SRC)))
INC := $(wildcard src/*.h)


test: test.exe
	@cd bin; echo "running program... \n Output Results:" ; mpirun ./test.exe;

Ctest: Ctest.exe
	@cd bin; echo 'running program...\nOutput Results:'; ./Ctest.exe;


lib: $(LIBDIR)/libPIC.a

$(LIBDIR)/libPIC.a: $(OBJ) 
	@echo make lib...
	ar rv $@ $^
	ranlib $@

%.exe: $(BINDIR)/%.o $(LIBDIR)/libPIC.a 
	@echo compilink and linking... 
	$(CC) -I src $< -o $(BINDIR)/$@ -L lib -l PIC

$(BINDIR)/%.o: %.C | $(INC)
	@echo compiling... $<
	$(CC) -I src -c $< -o $@


######### clean

tilde := $(wildcard */*~) $(wildcard *~)
exe := $(wildcard */*.exe) $(wildcard *.exe)
obj := $(wildcard */*.o) $(wildcard *.o)  $(wildcard */*.pcm) $(wildcard */*.d)
mylibs := $(wildcard */*.so) $(wildcard */*.a)

clean:
	@echo cleaning dir...
	rm -f $(exe) $(obj) $(tilde) $(mylibs)
