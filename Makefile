CC=g++
INC=-Iinclude
CFLAGS=-lncurses
DEPS=
LIB_GFX=src/cursesGfx.o src/cursesClock.o src/cursesGfx3d.o

all: cube chaosClocks testAsciiLines test3D stopwatch

%.o: %.cpp $(DEPS)
	$(CC) $(INC) -c -o $@ $<

stopwatch: src/stopwatch.o $(LIB_GFX)
	$(CC) $(input) $? -o $@  $(CFLAGS)

chaosClocks: src/chaosClocks.o $(LIB_GFX)
	$(CC) $(input) $? -o $@  $(CFLAGS)

testAsciiLines: src/testAsciiLines.o $(LIB_GFX)
	$(CC) $(input) $? -o $@  $(CFLAGS)

test3D: src/test3D.o $(LIB_GFX)
	$(CC) $(input) $? -o $@  $(CFLAGS)
	
cube: src/cube.o $(LIB_GFX)
	$(CC) $(input) $? -o $@  $(CFLAGS)
