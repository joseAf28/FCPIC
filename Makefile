# Makefile
BINDIR := bin
LIBDIR := lib

CCFLAGS :=  -Wextra -Wfloat-equal -Wundef -Werror -fverbose-asm  -Wshadow -Wpointer-arith -Wcast-align -Wconversion
DEBUG := -g -pedantic -Wall

CC :=  mpic++ -std=c++14 -O2

VPATH = main:src

SRC := $(wildcard src/*.cc)
OBJ := $(patsubst %.cc, $(BINDIR)/%.o, $(notdir $(SRC)))
INC := $(wildcard src/*.hh)

LIBH5 := -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial
INCH5 := -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -pthread

mpiAd: MPIAdvance.exe
	@cd bin; mpirun -np 3 ./MPIAdvance.exe -infile=test.txt

lib: $(LIBDIR)/libPIC.a

$(LIBDIR)/libPIC.a: $(OBJ) 
	@echo make lib...
	@ar rv $@ $^
	@ranlib $@

%.exe: $(BINDIR)/%.o $(LIBDIR)/libPIC.a 
	@echo compilink and linking... 
	@$(CC) -I src $< -o $(BINDIR)/$@ -L lib -l PIC $(LIBH5)


$(BINDIR)/%.o: %.cc | $(INC)
	@echo compiling... $<
	@$(CC) $(INCH5) -I src -c $< -o $@

######### clean

tilde := $(wildcard */*~) $(wildcard *~)
exe := $(wildcard */*.exe) $(wildcard *.exe)
obj := $(wildcard */*.o) $(wildcard *.o)  $(wildcard */*.pcm) $(wildcard */*.d)
mylibs := $(wildcard */*.so) $(wildcard */*.a)

txtOut := $(wildcard results/*/*.txt)
pngOut := $(wildcard results/*/*.png)
txtOutb := $(wildcard results/*.txt)
h5Out := $(wildcard results/*.h5)

clean:
	@echo cleaning dir...
	rm -f $(exe) $(obj) $(tilde) $(mylibs)

OutClean:
	@echo cleaning results...
	@rm -f $(txtOut) $(pngOut) $(txtOutb) $(h5Out)