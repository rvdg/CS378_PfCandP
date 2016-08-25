HOME       := /Users/rvdg
BLAS_LIB   := $(HOME)/blis/lib/libblis.a
CBLAS_INC  := $(HOME)/blis/include/blis
LAPACK_LIB := # $(HOME)/flame/lapack-3.3.0/lapack_LINUX.a
FLAME_LIB  := $(HOME)/flame/lib/libflame.a
FLAME_INC  := $(HOME)/flame/include

# indicate where the object files are to be created
CC         := gcc
FC         := gfortran
LINKER     := $(CC)
CFLAGS     := -O3 -Wall -I$(FLAME_INC) -I$(CBLAS_INC) -m64 -msse3
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed
NREPEATS   := 3#       number of times each experiment is repeated.
NFIRST     := 500#     smallest size to be timed.
NLAST      := 6000#    largest size to be timed.
NINC       := 500#     increment between sizes.

LDFLAGS    := -lm -lpthread -m64

# indicate where the FLAME include files reside

TEST_OBJS  := driver.o SymMatVec1.o SymMatVec2.o SymMatVec3.o

# $%.o: %.c
#	$(CC) $(CFLAGS) -c $< -o $@
# $%.o: %.f
#	$(FC) $(FFLAGS) -c $< -o $@

driver.x: $(TEST_OBJS)  $(LAPACK_OBJS)
	$(LINKER) $(TEST_OBJS) $(LAPACK_OBJS) $(LDFLAGS) $(FLAME_LIB) $(LAPACK_LIB) $(BLAS_LIB) -o driver.x

test:   driver.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver.x > output.m

clean:
	rm -f *.o *~ core *.x output.m


