# Makefile for twoDdipRY (2D Dipolar Hard Disk Integral Equation Solver)

# Compiler selection (default to ifx, can override via: make FC=gfortran)
FC = ifx

# Compiler-specific flags
ifeq ($(FC), ifx)
    FFLAGS = -O3 -qmkl
    LDFLAGS = -qmkl
else ifeq ($(FC), ifort)
    FFLAGS = -O3 -mkl
    LDFLAGS = -mkl
else ifeq ($(FC), gfortran)
    FFLAGS = -O3
    LDFLAGS = -lblas -llapack
endif

TARGET = 2DdipRYN_final
SRC = 2DdipRYN_final.f90

all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o *.mod

.PHONY: all clean
