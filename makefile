######################################################
# Makefile for sparsenc
# Ye Li
# leeyee.seu@gmail.com
######################################################

TOP = .
SRCDIR := .
OBJDIR := .
INCLUDEDIR = .
INC_PARMS = $(INCLUDEDIR:%=-I%)

UNAME := $(shell uname)
CC := gcc
ifeq ($(UNAME), Darwin)
	SED = gsed
	CC  = gcc-9
	#CC  = clang
	HAS_SSSE3 := $(shell sysctl -a | grep supplementalsse3)
	HAS_AVX2  := $(shell sysctl -a | grep avx2)
endif
ifeq ($(UNAME), Linux)
	SED = sed
	CC  = gcc
	HAS_SSSE3 := $(shell grep -i ssse3 /proc/cpuinfo)
	HAS_AVX2  := $(shell grep -i avx2 /proc/cpuinfo)
endif

CFLAGS0 = -Winline -std=c99 -lm -O3 -DNDEBUG $(INC_PARMS)
ifneq ($(HAS_SSSE3),)
	CFLAGS1 = -mssse3 -DINTEL_SSSE3
endif
ifneq ($(HAS_AVX2),)
	CFLAGS1 += -mavx2 -DINTEL_AVX2
endif
# Additional compile options
# CFLAGS2 = 

vpath %.h src include
vpath %.c src examples

DEFS    := galois.h bipartite.h bats.h bats-decoder-oa.h tiles.h channel.h
BATS-DYNBTS-SP    := $(OBJDIR)/galois.o $(OBJDIR)/bipartite.o $(OBJDIR)/bats-encoder.o $(OBJDIR)/bats-recoder.o $(OBJDIR)/mt19937ar.o $(OBJDIR)/gaussian.o $(OBJDIR)/bats-decoder-straight.c
$(OBJDIR)/%.o : $(OBJDIR)/%.c $(DEFS)
	$(CC) -c -o $@ $< $(CFLAGS0) $(CFLAGS1)
static-snc-Tp : $(BATS-DYNBTS-SP) static-bats-n-hop-Tp.c channel.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS1) $^ -lm
static-snc-Tp-fast : $(BATS-DYNBTS-SP) static-bats-n-hop-Tp-fast.c channel.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS1) $^ -lm
Q-learning-dynsnc-Tp : $(BATS-DYNBTS-SP) dynsnc-n-hop-Tp-Q-learning.c learning_functions.c channel.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS1) $^ -lm
MonteCarlo-dynsnc-Tp : $(BATS-DYNBTS-SP) dynsnc-n-hop-Tp-Monte-Carlo.c learning_functions.c channel.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS1) $^ -lm
Q-learning-dynsnc-Tp-fast : $(BATS-DYNBTS-SP) dynsnc-n-hop-Tp-Q-learning-fast.c learning_functions.c channel.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS1) $^ -lm

.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o Q-learning-dynsnc-Tp static-snc-Tp Q-learning-dynsnc-Tp-fast static-snc-Tp-fast MonteCarlo-dynsnc-Tp
