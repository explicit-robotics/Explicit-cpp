# [REF] https://stackoverflow.com/questions/30573481/how-to-write-a-makefile-with-separate-source-and-header-directories

SRC_DIR := src
OBJ_DIR := build
BIN_DIR := bin
# TEST_DIR := test

EXE  := $(BIN_DIR)/main
# TEST := $(BIN_DIR)/test
SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

CPPFLAGS := -Iinclude -I/opt/homebrew/include/eigen3 -MMD -MP
CFLAGS   := -Wall -g -std=c++17
CC = g++ 

.PHONY: all clean

all: $(EXE) #$(TEST)

$(EXE): $(OBJ) | $(BIN_DIR)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

# $(TEST): $(OBJ) | $(BIN_DIR)
# 	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv $(BIN_DIR) $(OBJ_DIR)

-include $(OBJ:.o=.d)