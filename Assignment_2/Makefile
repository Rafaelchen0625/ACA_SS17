
BIN=hevc_idct16
CC=gcc
CFLAGS=-O3 -msse4.2 -fno-tree-vectorize -std=gnu99
LIBS=-lrt

all: $(BIN)

hevc_idct16: hevc_idct16.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)


clean:
	rm $(BIN)
