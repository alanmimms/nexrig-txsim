# Makefile for NexRig Harmonic Analyzer

CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra
TARGET = harmonics
SOURCE = harmonics.cpp

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE) -lm

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) harmonic_analysis.csv

help:
	@echo "NexRig Harmonic Analyzer Build System"
	@echo ""
	@echo "Targets:"
	@echo "  make        - Build the analyzer"
	@echo "  make run    - Build and run the analyzer"
	@echo "  make clean  - Remove build artifacts"
	@echo "  make help   - Show this help message"
