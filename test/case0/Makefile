############################################
######### MAKEFILE FOR TEST CASE 0 #########
############################################

INCLUDE_DIR_GENERAL := ./../../include
INCLUDE_DIR_RBF := ./../../include/RBF

######################################################################################
############### THE USER NEED TO CHANGE HERE THE PATH TO EIGEN LIBRARY ###############
######################################################################################
INCLUDE_DIR_EIGEN := /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/include/eigen3
######################################################################################

SRC_DIR_RBF := ./../../src/RBF
OBJ_DIR := ./obj
BIN_DIR := ./bin

MAIN_SRC := ./main.cpp
RBF_FUNCTIONS_SRC := $(SRC_DIR_RBF)/RBFunctions.cpp

MAIN_OBJ := $(OBJ_DIR)/main.o
RBF_FUNCTIONS_OBJ := $(OBJ_DIR)/RBFunctions.o

EXECUTABLE := $(BIN_DIR)/prog

# Compilation
$(EXECUTABLE): $(MAIN_OBJ) $(RBF_FUNCTIONS_OBJ)
	@echo "Linking..."
	@mkdir -p $(BIN_DIR)
	@g++ -o $@ $^
	@echo "Done."	

$(OBJ_DIR)/main.o: $(MAIN_SRC)
	@echo "Compiling $<..."
	@mkdir -p $(OBJ_DIR)
	@g++ -o $@ -c $< -I$(INCLUDE_DIR_GENERAL) -I$(INCLUDE_DIR_RBF) -I$(INCLUDE_DIR_EIGEN)
	@echo "Done."

$(OBJ_DIR)/RBFunctions.o: $(RBF_FUNCTIONS_SRC)
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