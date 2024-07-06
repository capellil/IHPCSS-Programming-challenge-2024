CC=mpifort
CFLAGS=-O2 -Wall -mp=gpu -fopenmp

BIN_DIR=bin

default: all

all: create_directory \
	 $(BIN_DIR)/main

create_directory:
	@mkdir -p $(BIN_DIR)

bin/main: src/main.f08
	$(CC) -o $@ $< $(CFLAGS)
