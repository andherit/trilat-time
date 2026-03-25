FC = gfortran

# Keep flags modest by default; override from CLI if needed.
#FFLAGS ?= -O2 -Wall -Wextra
FFLAGS ?= -O3 -Wall -Wextra -march=native -ffast-math -unroll-loops

SRC_DIR := ../../src
OBJ_DIR := build/obj
MOD_DIR := build/mod

TARGET := two_layers

GENERIC_SRC  := $(SRC_DIR)/generic.f90
LISTS_SRC    := $(SRC_DIR)/lists.f90
LAT_MESH_SRC := $(SRC_DIR)/LAT_mesh.f90
LAT_TIME_SRC := $(SRC_DIR)/LAT_time.f90
TIME_SRC     := $(SRC_DIR)/time.f90
IOMOD_SRC    := $(SRC_DIR)/iomod.f90
MAIN_SRC     := two_layers.f90

GENERIC_OBJ  := $(OBJ_DIR)/generic.o
LISTS_OBJ    := $(OBJ_DIR)/lists.o
LAT_MESH_OBJ := $(OBJ_DIR)/LAT_mesh.o
LAT_TIME_OBJ := $(OBJ_DIR)/LAT_time.o
TIME_OBJ     := $(OBJ_DIR)/time.o
IOMOD_OBJ    := $(OBJ_DIR)/iomod.o
MAIN_OBJ     := $(OBJ_DIR)/two_layers.o

OBJS := \
	$(GENERIC_OBJ) \
	$(LISTS_OBJ) \
	$(LAT_MESH_OBJ) \
	$(LAT_TIME_OBJ) \
	$(TIME_OBJ) \
	$(IOMOD_OBJ) \
	$(MAIN_OBJ)

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

$(OBJ_DIR) $(MOD_DIR):
	mkdir -p $@

$(GENERIC_OBJ): $(GENERIC_SRC) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(LISTS_OBJ): $(LISTS_SRC) $(GENERIC_OBJ) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(LAT_MESH_OBJ): $(LAT_MESH_SRC) $(GENERIC_OBJ) $(LISTS_OBJ) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(LAT_TIME_OBJ): $(LAT_TIME_SRC) $(GENERIC_OBJ) $(LAT_MESH_OBJ) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(TIME_OBJ): $(TIME_SRC) $(GENERIC_OBJ) $(LISTS_OBJ) $(LAT_MESH_OBJ) $(LAT_TIME_OBJ) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(IOMOD_OBJ): $(IOMOD_SRC) $(GENERIC_OBJ) $(LAT_MESH_OBJ) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(MAIN_OBJ): $(MAIN_SRC) $(GENERIC_OBJ) $(LISTS_OBJ) $(LAT_MESH_OBJ) $(LAT_TIME_OBJ) $(TIME_OBJ) $(IOMOD_OBJ) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR) $(TARGET)
