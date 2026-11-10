CXX      ?= g++
CXXFLAGS := -std=c++17 -O3 -march=native -Wall -Wextra -pthread -MMD -MP
BUILD    := build
include books.mk

# One binary per approach: src/solvers/solve_<SOLVER>.cpp is linked with the shared main.
#   make run SOLVER=couzo ARGS="1 100"
SOLVER   ?= grid
BIN      := $(BUILD)/solve_$(SOLVER)

all: $(BIN)

$(BUILD)/solve_%: $(BUILD)/main.o $(BUILD)/solvers/solve_%.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Builds on the newest submissions/*.csv by default (keeps the better packing per n).
#   make run                       -> re-solve all puzzles 1..1000
#   make run ARGS="5 20"           -> re-solve only n = 5..20
#   make run ARGS="--base f.csv"   -> build on a specific file instead of the newest
#   make run ARGS="--fresh"        -> ignore existing submissions
run: $(BIN)
	./$(BIN) $(ARGS)

# Download/convert the open packings into books/ (see tools/fetch_books.py).
books:
	python3 tools/fetch_books.py

# Run every book solver one after another; each builds on the previous result, so the end state
# is the best packing per n across all books. A solver that improves nothing writes no file.
run-books: $(BUILD)/solve_grid
	@ls submissions/*.csv >/dev/null 2>&1 || ./$(BUILD)/solve_grid >/dev/null
	@for b in couzo $(filter-out couzo,$(BOOKS)); do \
	  echo "=== $$b"; $(MAKE) --no-print-directory run SOLVER=$$b ARGS="$(ARGS)" | grep -E "^base|^solved|^no change|^wrote|SCORE|skipped"; \
	done

# make validate [FILE=submissions/x.csv]    (defaults to the newest submission)
FILE ?= $(shell ls -t submissions/*.csv 2>/dev/null | head -1)
validate: $(BUILD)/solve_grid
	./$(BUILD)/solve_grid validate $(FILE)

clean:
	rm -rf $(BUILD)

-include $(wildcard $(BUILD)/*.d $(BUILD)/solvers/*.d)

.PHONY: all run books run-books validate clean
