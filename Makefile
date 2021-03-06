# Hey Emacs, this is a -*- makefile -*-

#---------------- Compiler Options C ----------------
CC = gcc

CFLAGS  = -m64 
CFLAGS += -std=c99
CFLAGS += -pedantic
CFLAGS += -Wall
CFLAGS += -Wshadow
CFLAGS += -Wpointer-arith
CFLAGS += -Wcast-qual
CFLAGS += -Wstrict-prototypes
CFLAGS += -Wmissing-prototypes
CFLAGS += -lm

#----------------------- Files -----------------------
FILES  = libsim.c 
FILES += physics/*.c 
FILES += math/*.c 
FILES += utils/*.c

#---------------------- Config -----------------------
BINDIR  = ./build/
LIBDIR  = ./build/lib/
TESTDIR = ./build/tests/


# Targets:
all: build lib

buildclean: clean all

ctesting:
	$(CC) test.c $(CFLAGS) -o test

build:
	mkdir -p $(BINDIR)
	$(CC) main.c $(FILES) $(CFLAGS) -o $(BINDIR)libsim

clean:
	rm -rf $(BINDIR)

test:
	rm -rf $(TESTDIR)
	mkdir -p $(TESTDIR)
	$(CC) tests/test.c tests/integrator.test.c tests/utils.test.c $(FILES) $(CFLAGS) -o $(TESTDIR)runtests
	$(TESTDIR)runtests

lib:
	mkdir $(LIBDIR)
	$(CC) $(CFLAGS) $(FILES)
	mv *.o $(LIBDIR)
	cd $(LIBDIR); ld -shared *.o -o libsim.so
	cd $(LIBDIR); rm -f *.o

.PHONY: build clean test
