# Makefile for Helmholtz EOS White Dwarf Generator
# ================================================

# Compiler and flags
FC = gfortran
FFLAGS = -O2 -Wall

# Source files
SOURCES = helmholtz_wrapper.f90 helmholtz_library.F90 generate_wd_eos.f90

# Object files
OBJECTS = $(SOURCES:.f90=.o)
OBJECTS := $(OBJECTS:.F90=.o)

# Target executable
TARGET = generate_wd_eos

# Default target
all: $(TARGET)

# Build the main executable
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^

# Pattern rules for object files
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

# Generate EOS table (runs the program)
eos: $(TARGET)
	./$(TARGET)
	cp wd_eos_helmholtz.dat ../wd_eos_helmholtz.dat
	@echo "Generated wd_eos_helmholtz.dat with relativistic enthalpy"
	@echo "Copied to parent directory for resampling"

# Convert to XNS format using Python resampling script
resample: eos
	cd .. && python3 resample_helmholtz.py
	@echo "Generated XNS-compatible resampled EOS table"

# Clean up generated files
clean:
	rm -f $(OBJECTS) $(TARGET) wd_eos_helmholtz.dat

# Clean everything including resampled files
clean-all: clean
	rm -f ../wd_eos_helmholtz_resampled.dat

# Help target
help:
	@echo "Available targets:"
	@echo "  all       - Build the EOS generator (default)"
	@echo "  eos       - Generate white dwarf EOS table"
	@echo "  resample  - Generate EOS and convert to XNS format"
	@echo "  clean     - Remove object files and executable"
	@echo "  clean-all - Remove all generated files including resampled EOS"
	@echo "  help      - Show this help message"

.PHONY: all eos resample clean clean-all help