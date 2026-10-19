# Compiler settings
CXX = g++
CXXFLAGS = -O3 -Wall -std=c++17 -Isrc

# Gurobi setup
# Users can override these by running: make GUROBI_HOME=/path/to/gurobi GUROBI_VER=110
GUROBI_HOME ?= /home/optlab/SharedSSD/Downloads/gurobi1203/linux64
GUROBI_VER ?= 120
GUROBI_INC = -I$(GUROBI_HOME)/include
GUROBI_LIB = -L$(GUROBI_HOME)/lib -Wl,-rpath,$(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi$(GUROBI_VER) -lm -lpthread

# Directories
SRCDIR = src
OBJDIR = obj
BINDIR = .

# Files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SOURCES))
EXECUTABLE = $(BINDIR)/ce_scheduler

# Default target
all: $(EXECUTABLE)

# Link the executable
$(EXECUTABLE): $(OBJECTS)
	@echo "Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(GUROBI_LIB)

# Compile source files to object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(OBJDIR)
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) $(GUROBI_INC) -c $< -o $@

# Clean target to remove compiled files
clean:
	@echo "Cleaning up..."
	rm -rf $(OBJDIR) $(EXECUTABLE)

.PHONY: all clean