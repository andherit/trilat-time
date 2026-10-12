FC = gfortran

# Keep flags modest by default; override from CLI if needed.
#FFLAGS ?= -O2 -Wall -Wextra
FFLAGS ?= -O3 -Wall -Wextra -march=native -ffast-math -unroll-loops
GMSH_LIB_DIR ?= $(HOME)/software/gmsh-4.15.2-Linux64-sdk/lib
LDFLAGS ?= -L$(GMSH_LIB_DIR) -Wl,-rpath,$(GMSH_LIB_DIR)
LDLIBS ?= -lgmsh

SRC_DIR := ../../../src
OBJ_DIR := build/obj
MOD_DIR := build/mod

TARGETS := mesh_gradient gradient

GENERIC_SRC  := $(SRC_DIR)/generic.f90
LISTS_SRC    := $(SRC_DIR)/lists.f90
LAT_MESH_SRC := $(SRC_DIR)/LAT_mesh.f90
LAT_TIME_SRC := $(SRC_DIR)/LAT_time.f90
TIME_DEBUG_SRC := $(SRC_DIR)/time_debug.f90
TIME_SRC     := $(SRC_DIR)/time.f90
IOMOD_SRC    := $(SRC_DIR)/iomod.f90
GMSH_SRC     := gmsh.f90
MESH_MAIN_SRC     := mesh_gradient.f90
GRADIENT_MAIN_SRC := gradient.f90

GENERIC_OBJ  := $(OBJ_DIR)/generic.o
LISTS_OBJ    := $(OBJ_DIR)/lists.o
LAT_MESH_OBJ := $(OBJ_DIR)/LAT_mesh.o
LAT_TIME_OBJ := $(OBJ_DIR)/LAT_time.o
TIME_DEBUG_OBJ := $(OBJ_DIR)/time_debug.o
TIME_OBJ     := $(OBJ_DIR)/time.o
IOMOD_OBJ    := $(OBJ_DIR)/iomod.o
GMSH_OBJ     := $(OBJ_DIR)/gmsh.o
MESH_MAIN_OBJ     := $(OBJ_DIR)/mesh_gradient.o
GRADIENT_MAIN_OBJ := $(OBJ_DIR)/gradient.o

MESH_OBJS := \
	$(GENERIC_OBJ) \
	$(LISTS_OBJ) \
	$(LAT_MESH_OBJ) \
	$(IOMOD_OBJ) \
	$(GMSH_OBJ) \
	$(MESH_MAIN_OBJ)

GRADIENT_OBJS := \
	$(GENERIC_OBJ) \
	$(LISTS_OBJ) \
	$(LAT_MESH_OBJ) \
	$(LAT_TIME_OBJ) \
	$(TIME_DEBUG_OBJ) \
	$(TIME_OBJ) \
	$(IOMOD_OBJ) \
	$(GRADIENT_MAIN_OBJ)

.PHONY: all run clean

all: $(TARGETS)

mesh_gradient: $(MESH_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(MESH_OBJS) $(LDLIBS)

gradient: $(GRADIENT_OBJS)
	$(FC) $(FFLAGS) -o $@ $(GRADIENT_OBJS)

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

$(TIME_DEBUG_OBJ): $(TIME_DEBUG_SRC) $(GENERIC_OBJ) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(TIME_OBJ): $(TIME_SRC) $(GENERIC_OBJ) $(LISTS_OBJ) $(LAT_MESH_OBJ) $(LAT_TIME_OBJ) $(TIME_DEBUG_OBJ) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(IOMOD_OBJ): $(IOMOD_SRC) $(GENERIC_OBJ) $(LAT_MESH_OBJ) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(GMSH_OBJ): $(GMSH_SRC) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(MESH_MAIN_OBJ): $(MESH_MAIN_SRC) $(GENERIC_OBJ) $(LISTS_OBJ) $(LAT_MESH_OBJ) $(IOMOD_OBJ) $(GMSH_OBJ) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(GRADIENT_MAIN_OBJ): $(GRADIENT_MAIN_SRC) $(GENERIC_OBJ) $(LISTS_OBJ) $(LAT_MESH_OBJ) $(LAT_TIME_OBJ) $(TIME_DEBUG_OBJ) $(TIME_OBJ) $(IOMOD_OBJ) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

run: gradient
	./gradient

clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR) $(TARGETS)
