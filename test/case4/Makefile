# Variables
INCLUDE_DIR_GENERAL := ./../../include
INCLUDE_DIR_RBF := ./../../include/RBF
INCLUDE_DIR_OLS := ./../../include/OLS


################################################################
#### THE USER NEED TO CHANGE HERE THE PATH TO EIGEN LIBRARY ####
################################################################
INCLUDE_DIR_EIGEN := /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/include/eigen3

SRC_DIR_RBF := ./../../src/RBF
SRC_DIR_OLS := ./../../src/OLS
SRC_DIR_RESCALE := ./../../src
OBJ_DIR := ./obj
BIN_DIR := ./bin

MAIN_SRC := ./main.cpp
RBF_INTERPOLATOR_SRC := $(SRC_DIR_RBF)/RBFinterpolator.cpp
RBF_FUNCTIONS_SRC := $(SRC_DIR_RBF)/RBFunctions.cpp
OLS_INTERPOLATOR_SRC := $(SRC_DIR_OLS)/OLSinterpolator.cpp
RESCALE_SRC := $(SRC_DIR_RESCALE)/rescaling.cpp

MAIN_OBJ := $(OBJ_DIR)/main.o
RBF_INTERPOLATOR_OBJ := $(OBJ_DIR)/RBFinterpolator.o
RBF_FUNCTIONS_OBJ := $(OBJ_DIR)/RBFunctions.o
OLS_INTERPOLATOR_OBJ := $(OBJ_DIR)/OLSinterpolator.o
RESCALE_OBJ := $(OBJ_DIR)/rescaling.o

EXECUTABLE := $(BIN_DIR)/prog

# Compilation
$(EXECUTABLE): $(MAIN_OBJ) $(RBF_INTERPOLATOR_OBJ) $(RBF_FUNCTIONS_OBJ) $(OLS_INTERPOLATOR_OBJ) $(RESCALE_OBJ)
	@echo "Linking..."
	@mkdir -p $(BIN_DIR)
	@g++ -o $@ $^
	@echo "Done."	

$(OBJ_DIR)/main.o: $(MAIN_SRC)
	@echo "Compiling $<..."
	@mkdir -p $(OBJ_DIR)
	@g++ -o $@ -c $< -I$(INCLUDE_DIR_GENERAL) -I$(INCLUDE_DIR_RBF) -I$(INCLUDE_DIR_OLS) -I$(INCLUDE_DIR_EIGEN)
	@echo "Done."

$(OBJ_DIR)/RBFinterpolator.o: $(RBF_INTERPOLATOR_SRC)
	@echo "Compiling $<..."
	@mkdir -p $(OBJ_DIR)
	@g++ -o $@ -c $< -I$(INCLUDE_DIR_GENERAL) -I$(INCLUDE_DIR_RBF) -I$(INCLUDE_DIR_OLS) -I$(INCLUDE_DIR_EIGEN)
	@echo "Done."

$(OBJ_DIR)/RBFunctions.o: $(RBF_FUNCTIONS_SRC)
	@echo "Compiling $<..."
	@mkdir -p $(OBJ_DIR)
	@g++ -o $@ -c $< -I$(INCLUDE_DIR_GENERAL) -I$(INCLUDE_DIR_RBF) -I$(INCLUDE_DIR_OLS) -I$(INCLUDE_DIR_EIGEN)
	@echo "Done."

$(OBJ_DIR)/OLSinterpolator.o: $(OLS_INTERPOLATOR_SRC)
	@echo "Compiling $<..."
	@mkdir -p $(OBJ_DIR)
	@g++ -o $@ -c $< -I$(INCLUDE_DIR_GENERAL) -I$(INCLUDE_DIR_RBF) -I$(INCLUDE_DIR_OLS) -I$(INCLUDE_DIR_EIGEN)
	@echo "Done."

$(OBJ_DIR)/rescaling.o: $(RESCALE_SRC)
	@echo "Compiling $<..."
	@mkdir -p $(OBJ_DIR)
	@g++ -o $@ -c $< -I$(INCLUDE_DIR_GENERAL) -I$(INCLUDE_DIR_RBF) -I$(INCLUDE_DIR_EIGEN)
	@echo "Done."

# run
run:
	@./bin/prog


# Clean
clean:
	@echo "Cleaning..."
	@rm -rf $(OBJ_DIR) $(BIN_DIR) ./plot/
	@echo "Done."

# Phony targets
.PHONY: clean