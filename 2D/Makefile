# Compiles by default with the DEBUG flag
CFLAGS = -DDEBUG

all: ./main.o
clean:
	find . -iname '*.o' -print -delete

./main.o: main.cpp ./include/utils.o
	g++ $(CFLAGS) main.cpp ./include/utils.o -o ./main.o
./include/utils.o: include/utils.h include/utils.c
	g++ $(CFLAGS) -c include/utils.c -o ./include/utils.o
