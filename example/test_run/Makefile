PROGNAME := nicole

# Define source directories
SRC_DIR := ../../nicole_src
UTILS_DIR := utils
ELEMENT_DIR := elements
SPECIES_DIR := species
REACTION_DIR := reactions
SIMULATOR_DIR := reaction_simulator
RESISTIVITY_DIR := mhd_resistivity
ODEPACK_DIR := odepack_cpp

# Create build directory under test_run
OUT_DIR := build

# Define the final target (executable)
TARGET := $(OUT_DIR)/$(PROGNAME)

# Collect all source files (relative paths to correct files)
SRCS := $(wildcard *.cpp) \
        $(wildcard $(SRC_DIR)/*.cpp) \
        $(wildcard $(SRC_DIR)/$(UTILS_DIR)/*.cpp) \
        $(wildcard $(SRC_DIR)/$(ELEMENT_DIR)/*.cpp) \
        $(wildcard $(SRC_DIR)/$(SPECIES_DIR)/*.cpp) \
        $(wildcard $(SRC_DIR)/$(REACTION_DIR)/*.cpp) \
        $(wildcard $(SRC_DIR)/$(SIMULATOR_DIR)/*.cpp) \
        $(wildcard $(SRC_DIR)/$(RESISTIVITY_DIR)/*.cpp) \
        $(wildcard $(SRC_DIR)/$(ODEPACK_DIR)/*.cpp)

# Set the paths for the object files in the build directory
OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(OUT_DIR)/%.o)

# Define the compiler
CPP = g++
# Define compilation flags
CFLAGS = -g -c -O3 -std=c++17

# Declare all and clean as phony targets
.PHONY: all clean

# Default target to build the executable
all: $(TARGET)

# Linking rule for creating the executable from object files
$(TARGET): $(OBJS)
	$(CPP) -o $@ $^

# Rule for compiling source files into object files
$(OUT_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)  # Create necessary subdirectories in build directory
	$(CPP) $(CFLAGS) -o $@ -c $<

# Clean up the build directory and executable
clean:
	rm -rf $(OUT_DIR) $(TARGET)
