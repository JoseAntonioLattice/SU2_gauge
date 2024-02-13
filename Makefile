FC = gfortran

SRC=src
BIN=bin

TARGET = SU2_gauge.exe

SOURCE = parameters.f90 data_types_observables.f90 arrays.f90 starts.f90 dynamics.f90 main.f90

OBJECT = $(patsubst %,$(BIN)/%, $(notdir  $(SOURCE:.f90=.o)))

FFLAGS = -J$(BIN) -I$(BIN) -fcheck=all -fbacktrace -O0 -Wall -Wextra -std=f2008

$(BIN)/$(TARGET): $(OBJECT)
	$(FC) -o $@ $^

$(BIN)/%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY: help run clean

run:
	echo "input_parameters.par" | $(BIN)/$(TARGET)

help:
	@echo "src: $(SOURCE)"
	@echo "bin: $(OBJECT)"

clean :
	rm -f $(OBJECT) $(BIN)/$(TARGET)
